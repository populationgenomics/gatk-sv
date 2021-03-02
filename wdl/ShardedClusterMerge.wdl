version 1.0

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks

# Merges
workflow ShardedClusterMerge {
  input {
    Array[File] vcfs
    String id_prefix
    String file_prefix

    String sv_pipeline_docker
    String sv_base_mini_docker

    # Do not use
    File? NONE_FILE_

    # overrides for local tasks
    RuntimeAttr? runtime_override_merge_clusters

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_concat_shards
  }

  #Run vcfcluster per shard
  scatter (i in range(length(vcfs))) {
    call MergeClusters {
      input:
        vcf=vcfs[i],
        shard=i,
        file_prefix=file_prefix,
        id_prefix=id_prefix,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_merge_clusters
    }
  }

  call MiniTasks.ConcatVcfs as ConcatShards {
    input:
      vcfs=MergeClusters.out,
      vcfs_idx=MergeClusters.out_index,
      merge_sort=true,
      outfile_prefix=file_prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_shards
  }

  #Output
  output {
    File clustered_vcf = ConcatShards.concat_vcf
    File clustered_vcf_idx = ConcatShards.concat_vcf_idx
  }
}

#Run svtk vcfcluster
task MergeClusters {
  input {
    File vcf
    String file_prefix
    String id_prefix
    Int shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_name = file_prefix + ".merge_shard_" + shard

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float shard_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0,
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
    memory: runtime_mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    gunzip -c ~{vcf} > unmerged.vcf
    echo "unmerged.vcf" > vcf.list

    svtk vcfcluster vcf.list merged.vcf \
      -p ~{id_prefix} \
      --preserve-ids \
      --preserve-genotypes \
      --preserve-header \
      --merge-only

    # remove CLUSTER field, sort and compress the vcf
    mkdir temp
    bcftools annotate --no-version -x INFO/CLUSTER merged.vcf \
      | bcftools sort --temp-dir temp --max-mem ~{sort_mem_mb} --output-type z --output-file ~{output_name}
    tabix ~{output_name}
  >>>

  output {
    File clustered_vcf = output_name
    File clustered_vcf_idx = output_name + ".tbi"
  }
}
