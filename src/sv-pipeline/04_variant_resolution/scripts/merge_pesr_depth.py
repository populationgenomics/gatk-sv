#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
from collections import defaultdict
import pysam
import svtk.utils as svu


def merge_pesr_depth(vcf, fout, prefix, frac=0.5, sample_overlap=0.5):

    # Given one pesr record and one depth record, merge depth attributes into the pesr record
    def _merge_pair(record_a, record_b):
        is_depth_a = _record_is_depth(record_a)
        is_depth_b = _record_is_depth(record_b)
        if is_depth_a == is_depth_b:
            raise ValueError("Attempted to write pesr/pesr or depth/depth pair")
        if is_depth_a:
            depth_record = record_a
            pesr_record = record_b
        else:
            pesr_record = record_a
            depth_record = record_b

        pesr_record.info['ALGORITHMS'] = tuple(sorted(list(set(pesr_record.info['ALGORITHMS'] + ('depth',)))))
        pesr_record.info['MEMBERS'] = (pesr_record.id, depth_record.id)

        for sample in pesr_record.samples:
            if 'EV' in pesr_record.samples[sample].keys() and 'EV' in depth_record.info.keys():
                pesr_ev = pesr_record.samples[sample]['EV']
                depth_ev = depth_record.samples[sample]['EV']
                pesr_record.samples[sample]['EV'] = tuple(sorted(set(pesr_ev).union(depth_ev)))

        if 'varGQ' in pesr_record.info.keys() and 'varGQ' in depth_record.info.keys():
            pesr_record.info['varGQ'] = max(pesr_record.info['varGQ'], depth_record.info['varGQ'])

    def _write_record(record, salvaged):
        svtype = record.info['SVTYPE']
        counts[svtype] += 1
        if _record_is_depth(record):
            new_record = clean_depth_record(base_record, record)
        else:
            new_record = record.copy()
        member_id = new_record.id
        new_record.id = '{0}_{1}_{2}'.format(prefix, svtype, counts[svtype])
        if salvaged:
            new_record.id = new_record.id + "_salvaged"
        new_record.info['MEMBERS'] = (member_id,)
        fout.write(new_record)

    def _flush_active_records():
        for r in active_records:
            _write_record(r, False)
        active_records.clear()

    def _record_is_depth(record):
        alg = record.info['ALGORITHMS']
        return len(alg) == 1 and alg[0] == 'depth'

    def _reciprocal_overlap(record_a, record_b):
        return svu.recip(record_a.pos, record_a.stop, record_b.pos, record_b.stop, frac=frac)

    def _sample_overlap(record_a, record_b):
        return svu.samples_overlap(svu.get_called_samples(record_a), svu.get_called_samples(record_b), upper_thresh=sample_overlap, lower_thresh=sample_overlap)

    def _records_cluster_together(record_a, record_b):
        return _record_is_depth(record_a) != _record_is_depth(record_b) \
               and record_a.info['SVTYPE'] == record_b.info['SVTYPE'] \
               and _reciprocal_overlap(record_a, record_b) \
               and _sample_overlap(record_a, record_b)

    def _get_base_record(vcf):
        for record in vcf.fetch():
            if not _record_is_depth(record):
                vcf.reset()
                return record

    base_record = _get_base_record(vcf)
    if base_record is None:
        raise ValueError("No PESR records were found")

    active_records = []
    counts = defaultdict(int)
    current_contig = None

    for record in vcf.fetch():

        # Write all-ref sites as "salvaged"
        samples = svu.get_called_samples(record)
        if len(samples) == 0:
            _write_record(record, True)
            continue

        finalized_record_ids = set()
        clustered_depth_ids = set()
        if current_contig is None or record.contig != current_contig:
            # Started a new contig
            _flush_active_records()
            current_contig = record.chrom
        else:
            # Check if this record belongs to an existing cluster
            for ar in active_records:
                if ar.stop < record.start:
                    # Since traversing in order, this cluster cannot have any more members
                    # Only write depth record if not clustered with a pesr record (filtered in comprehension below)
                    _write_record(ar, False)
                    finalized_record_ids.add(ar.id)
                elif _records_cluster_together(record, ar):
                    _merge_pair(record, ar)
                    if _record_is_depth(record):
                        clustered_depth_ids.add(record.id)
                    if _record_is_depth(ar):
                        clustered_depth_ids.add(ar.id)
        active_records.append(record)
        active_records = [r for r in active_records if r.id not in finalized_record_ids and r.id not in clustered_depth_ids]

    _flush_active_records()


def clean_depth_record(base_record, depth_record):
    base = base_record.copy()
    base.chrom = depth_record.chrom
    base.pos = depth_record.pos
    base.id = depth_record.id
    base.ref = depth_record.ref
    base.alts = depth_record.alts
    base.stop = depth_record.stop

    for key in base.info.keys():
        if key not in depth_record.info:
            base.info.pop(key)

    for key, val in depth_record.info.items():
        base.info[key] = val

    for sample in depth_record.samples:
        for key, val in depth_record.samples[sample].items():
            base.samples[sample][key] = val

    return base


def check_header(vcf):
    if 'MEMBERS' not in vcf.header.info:
        vcf.header.add_line("##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Combined but unmerged VCF of PE/SR calls')
    parser.add_argument('fout', help='Output VCF (unsorted!), can be "-" or "stdout"')
    parser.add_argument('--prefix', default='pesr_rd_merged')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    check_header(vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    merge_pesr_depth(vcf, fout, args.prefix)
    fout.close()


if __name__ == '__main__':
    main()
