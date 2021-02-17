version 1.0

import "Structs.wdl"
import "GATKSVPreprocessSample.wdl" as preprocess

workflow GATKSVPreprocessSampleScatter {
  input {
    Array[String] sample_id
    Array[File] manta_vcf
    Array[File] melt_vcf
    Array[File] wham_vcf
    Array[Array[File]] cnv_beds
    # TODO : gcnv vcfs
  }

  scatter (i in range(length(sample_id))) {
    call preprocess.GATKSVPreprocessSample {
      input:
        sample_id = sample_id[i],
        manta_vcf = manta_vcf[i],
        melt_vcf = melt_vcf[i],
        wham_vcf = wham_vcf[i],
        cnv_beds = cnv_beds[i]
    }
  }

  output {
    Array[File] out = GATKSVPreprocessSample.out
    Array[File] out_index = GATKSVPreprocessSample.out_index
  }
}