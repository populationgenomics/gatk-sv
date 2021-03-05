version 1.0
# based on snapshot 6
# https://portal.firecloud.org/#methods/Talkowski-SV/04_sharded_vcfcluster/6/wdl

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks
import "ShardedClusterMerge.wdl" as scm

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
    Float merging_shard_scale_factor = 500000000

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
  }

  File vcf_idx = vcf + ".tbi"
  if (defined(exclude_list)) {
    File exclude_list_idx = exclude_list + ".tbi"
  }
  String sv_type_prefix = prefix + "." + contig + "." + sv_type

  #New as of November 2, 2018: perform sharding and return list of variant IDs
  # for each shard, rather than VCF shards themselves, which should dramatically
  # improve speed of sharding task (previously took 1-6 hours for 14k samples in
  # gnomAD v2)
  call ShardVcfPrecluster {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      dist=dist,
      frac=frac,
      max_shards=max_shards,
      min_per_shard=min_per_shard,
      prefix=sv_type_prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_shard_vcf_precluster
  }

  Int merge_shard_size = ceil(merging_shard_scale_factor / ShardVcfPrecluster.num_samples)

  #Run vcfcluster per shard
  scatter (i in range(length(ShardVcfPrecluster.VID_list_shards))) {
    call SvtkVcfCluster {
      input:
        vcf=vcf,
        shard=i,
        VIDs=ShardVcfPrecluster.VID_list_shards[i],
        prefix=sv_type_prefix,
        dist=dist,
        frac=frac,
        sample_overlap=sample_overlap,
        exclude_list=exclude_list,
        exclude_list_idx=exclude_list_idx,
        svsize=sv_size,
        records_per_shard=merge_shard_size,
        sv_types=sv_types,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_svtk_vcf_cluster
    }
    call scm.ShardedClusterMerge {
      input:
        vcfs = SvtkVcfCluster.out,
        id_prefix = sv_type_prefix,
        file_prefix = sv_type_prefix + ".vid_shard_" + i,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_override_merge_clusters = runtime_override_merge_clusters,
        runtime_override_concat_shards = runtime_override_concat_inner_shards,
        runtime_override_sort_merged_vcf = runtime_override_sort_merged_vcf
    }
  }
  if (length(ShardedClusterMerge.clustered_vcf) == 0) {
    call GetVcfHeaderWithMembersInfoLine {
      input:
        vcf_gz=vcf,
        prefix=sv_type_prefix,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_get_vcf_header_with_members_info_line
    }
  }

  #Merge shards per svtype
  Array[File] clustered_vcfs = if length(ShardedClusterMerge.clustered_vcf) > 0 then ShardedClusterMerge.clustered_vcf else select_all([GetVcfHeaderWithMembersInfoLine.out])
  Array[File] clustered_vcfs_idx = if length(ShardedClusterMerge.clustered_vcf) > 0 then ShardedClusterMerge.clustered_vcf_idx else select_all([GetVcfHeaderWithMembersInfoLine.out_idx])
  call MiniTasks.ConcatVcfs as ConcatShards {
    input:
      vcfs=clustered_vcfs,
      vcfs_idx=clustered_vcfs_idx,
      merge_sort=true,
      outfile_prefix=sv_type_prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_shards
  }

  #Output
  output {
    File clustered_vcf = ConcatShards.concat_vcf
    File clustered_vcf_idx = ConcatShards.concat_vcf_idx
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

#Intelligently shard a VCF for parallelized clustering
task ShardVcfPrecluster {
  input {
    File vcf
    File vcf_idx
    Int dist
    Float frac
    Int max_shards
    Int min_per_shard
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float shard_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(base_disk_gb + shard_size * input_disk_scale),
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
    set -eu -o pipefail
    bcftools query -l ~{vcf} | wc -l > sample_count.txt

    /opt/sv-pipeline/04_variant_resolution/scripts/shardVCF_preClustering_part1.sh \
      -D ~{dist} \
      -R ~{frac} \
      -L ~{min_per_shard} \
      -S ~{max_shards} \
      -P ~{prefix} \
      ~{vcf}
  >>>

  output {
    Array[File] VID_list_shards = glob("*.VIDs.list")
    Int num_samples = read_int("sample_count.txt")
  }
}


#Run svtk vcfcluster
task SvtkVcfCluster {
  input {
    File vcf
    Int shard
    File VIDs
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

  String output_prefix = "~{prefix}.shard_~{shard}.unmerged_clusters"

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float shard_size = size(select_all([vcf, VIDs, exclude_list, exclude_list_idx]), "GiB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + shard_size * input_mem_scale,
    disk_gb: ceil(base_disk_gb + shard_size * input_disk_scale),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  Float runtime_mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
  Int sort_mem_mb = floor(runtime_mem_gb * 1000 - 100)
  runtime {
    memory: "~{runtime_mem_gb} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    bcftools view --include ID=@~{VIDs} ~{vcf} -O z -o input.vcf.gz
    tabix input.vcf.gz

    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -f -p bed ~{exclude_list}" else ""}

    #Run clustering
    echo "input.vcf.gz" > unclustered_vcfs.list
    svtk vcfcluster unclustered_vcfs.list unmerged_clusters.vcf \
      -d ~{dist} \
      -f ~{frac} \
      ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
      -z ~{svsize} \
      -p ~{prefix} \
      -t ~{sep=',' sv_types} \
      -o ~{sample_overlap} \
      --preserve-ids \
      --preserve-genotypes \
      --preserve-header \
      --skip-merge

    # Shard output
    python3 <<CODE
    import sys
    import gzip
    import os

    header = []
    with open('unmerged_clusters.vcf') as f:
      current_vid = ""
      current_cluster = []
      current_shard = 0
      current_shard_size = 0
      shard_path_format = "~{output_prefix}.cluster_shard_{}.vcf.gz"
      shard_path = shard_path_format.format(current_shard)
      fout = gzip.open(shard_path, 'wb')
      if fout is None:
        raise IOError("Could not open '{}'".format(shard_path))
        sys.exit(1)
      for line in f:
        if line[0] == '#':
          fout.write(line.encode('utf-8'))
          header.append(line)
        else:
          tok = line.split('\t', 3)
          vid = tok[2]
          if vid == current_vid:
            current_cluster.append(line)
          else:
            for record in current_cluster:
              fout.write(record.encode('utf-8'))
            current_shard_size += len(current_cluster)
            if current_shard_size >= ~{records_per_shard}:
              current_shard += 1
              current_shard_size = 0
              fout.close()
              shard_path = shard_path_format.format(current_shard)
              fout = gzip.open(shard_path, 'wb')
              if fout is None:
                raise IOError("Could not open '{}'".format(shard_path))
                sys.exit(1)
              for hline in header:
                fout.write(hline.encode('utf-8'))
            current_cluster = [line]
            current_vid = vid

      # Write last cluster
      for record in current_cluster:
        fout.write(record.encode('utf-8'))
      current_shard_size += len(current_cluster)
      fout.close()

      # Delete trailing empty shard
      if current_shard > 0 and current_shard_size == 0:
        os.remove(shard_path)
    CODE

  >>>

  output {
    # NOT block-compressed
    Array[File] out = glob("~{output_prefix}.cluster_shard_*.vcf.gz")
  }
}
