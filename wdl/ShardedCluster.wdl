version 1.0
# based on snapshot 6
# https://portal.firecloud.org/#methods/Talkowski-SV/04_sharded_vcfcluster/6/wdl

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks
import "Utils.wdl" as utils

# Workflow to shard a filtered vcf & run vcfcluster (sub-sub-sub workflow)
workflow ShardedCluster {
  input {
    File vcf
    Int dist
    Float frac
    Int max_shards
    Int min_per_shard
    String prefix
    String contig
    String sv_type
    Float sample_overlap
    File? exclude_list
    Int sv_size
    Array[String] sv_types
    Float merging_shard_scale_factor = 100000000

    String sv_pipeline_docker
    String sv_base_mini_docker

    # Do not use
    File? NONE_FILE_

    # overrides for local tasks
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line

    # overrides for merge subworkflow
    RuntimeAttr? runtime_override_merge_clusters
    RuntimeAttr? runtime_override_concat_inner_shards

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_concat_shards
    RuntimeAttr? runtime_override_sort_merged_vcf
    RuntimeAttr? runtime_override_count_samples
  }

  File vcf_idx = vcf + ".tbi"
  if (defined(exclude_list)) {
    File exclude_list_idx = exclude_list + ".tbi"
  }
  String sv_type_prefix = prefix + "." + contig + "." + sv_type

  call utils.CountSamples {
    input:
    vcf=vcf,
    sv_base_mini_docker=sv_base_mini_docker,
    runtime_attr_override=runtime_override_count_samples
  }
  Int merge_shard_size = ceil(merging_shard_scale_factor / CountSamples.num_samples)

  call ShardClusters {
    input:
      vcf=vcf,
      prefix=sv_type_prefix,
      dist=dist,
      frac=frac,
      exclude_list=exclude_list,
      exclude_list_idx=exclude_list_idx,
      svsize=sv_size,
      records_per_shard=merge_shard_size,
      sv_types=sv_types,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_svtk_vcf_cluster
  }

  #Run vcfcluster per shard
  scatter (i in range(length(ShardClusters.out))) {
    call SvtkVcfCluster {
      input:
        vcf=vcf,
        vids=ShardClusters.out[i],
        prefix="~{sv_type_prefix}.shard_${i}.clustered",
        dist=dist,
        frac=frac,
        exclude_list=exclude_list,
        exclude_list_idx=exclude_list_idx,
        svsize=sv_size,
        sample_overlap=sample_overlap,
        records_per_shard=merge_shard_size,
        sv_types=sv_types,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_svtk_vcf_cluster
    }
    call MiniTasks.SortVcf {
      input:
        vcf = SvtkVcfCluster.out,
        outfile_prefix = "~{sv_type_prefix}.shard_${i}.sorted",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_sort_merged_vcf
    }
  }

  if (length(SvtkVcfCluster.out) == 0) {
    call GetVcfHeaderWithMembersInfoLine {
      input:
        vcf_gz=vcf,
        prefix=sv_type_prefix,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_get_vcf_header_with_members_info_line
    }
  }
  if (length(SvtkVcfCluster.out) > 0) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=SortVcf.out,
        vcfs_idx=SortVcf.out_index,
        merge_sort=true,
        outfile_prefix="~{sv_type_prefix}.",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_shards
    }
  }

  #Output
  output {
    File clustered_vcf = select_first([GetVcfHeaderWithMembersInfoLine.out, ConcatVcfs.concat_vcf])
    File clustered_vcf_idx = select_first([GetVcfHeaderWithMembersInfoLine.out_idx, ConcatVcfs.concat_vcf_idx])
  }
}

# Adds MEMBERS definition to header (workaround for when VIDs_list is empty)
task GetVcfHeaderWithMembersInfoLine {
  input {
    File vcf_gz
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    gunzip -c ~{vcf_gz} | grep "^##" > header
    echo "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">" >> header
    gunzip -c ~{vcf_gz} | grep "^#" | grep -v "^##" >> header
    bgzip -c header > ~{prefix}.members.vcf.gz
    tabix ~{prefix}.members.vcf.gz
  >>>

  output {
    File out = "~{prefix}.members.vcf.gz"
    File out_idx = "~{prefix}.members.vcf.gz.tbi"
  }
}

#Do fast cluster on sites-only vcf (sample_overlap = 0) to generate shards
task ShardClusters {
  input {
    File vcf
    String prefix
    Int dist
    Float frac
    File? exclude_list
    File? exclude_list_idx
    Int svsize
    Array[String] sv_types
    Int records_per_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_prefix = "~{prefix}.unmerged_clusters"

  Float input_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 1.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    bcftools view -G ~{vcf} -o sites_only.vcf
    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}
    #Run clustering
    svtk vcfcluster <(echo "sites_only.vcf") unmerged_clusters.vcf \
      -d ~{dist} \
      -f ~{frac} \
      ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
      -z ~{svsize} \
      -p ~{prefix} \
      -t ~{sep=',' sv_types} \
      -o 0 \
      --preserve-header \
      --preserve-ids \
      --skip-merge

    # Shard output
    python3 <<CODE
    import sys
    import os

    with open('unmerged_clusters.vcf') as f:
      current_vid = ""
      current_cluster = []
      current_shard = 0
      current_shard_size = 0
      shard_path_format = "~{output_prefix}.vids.shard_{}.list"
      shard_path = shard_path_format.format(current_shard)
      fout = open(shard_path, 'w')
      if fout is None:
        raise IOError("Could not open '{}'".format(shard_path))
        sys.exit(1)
      for line in f:
        if line[0] == '#':
        else:
          tok = line.split('\t', 3)
          vid = tok[2]
          if vid == current_vid:
            current_cluster.append(line)
          else:
            for record in current_cluster:
              fout.write(record.id + '\n')
            current_shard_size += len(current_cluster)
            if current_shard_size >= ~{records_per_shard}:
              current_shard += 1
              current_shard_size = 0
              fout.close()
              shard_path = shard_path_format.format(current_shard)
              fout = open(shard_path, 'w')
              if fout is None:
                raise IOError("Could not open '{}'".format(shard_path))
                sys.exit(1)
            current_cluster = [line]
            current_vid = vid

      # Write last cluster
      for record in current_cluster:
        fout.write(record.id + '\n')
      current_shard_size += len(current_cluster)
      fout.close()

      # Delete trailing empty shard
      if current_shard > 0 and current_shard_size == 0:
        os.remove(shard_path)
    CODE

  >>>

  output {
    Array[File] out = glob("~{output_prefix}.vids.shard_*.list")
  }
}

task SvtkVcfCluster {
  input {
    File vcf
    File vids
    String prefix
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    File? exclude_list_idx
    Int svsize
    Array[String] sv_types
    Int records_per_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_prefix = "~{prefix}"

  Float input_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    bcftools view --no-version --include ID=@~{vids} ~{vcf} -O z -o unclustered.vcf.gz
    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}
    #Run clustering
    svtk vcfcluster <(echo "unclustered.vcf.gz") ~{output_prefix}.vcf \
      -d ~{dist} \
      -f ~{frac} \
      ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
      -z ~{svsize} \
      -p ~{prefix} \
      -t ~{sep=',' sv_types} \
      -o ~{sample_overlap} \
      --preserve-ids \
      --preserve-genotypes \
      --preserve-header
    gzip ~{output_prefix}.vcf
  >>>

  output {
    File out = "~{output_prefix}.vcf.gz"
  }
}