"use strict";(self.webpackChunkGATK_SV=self.webpackChunkGATK_SV||[]).push([[5684],{3447:(e,i,l)=>{l.r(i),l.d(i,{assets:()=>d,contentTitle:()=>a,default:()=>p,frontMatter:()=>o,metadata:()=>n,toc:()=>c});const n=JSON.parse('{"id":"modules/filter_genotypes","title":"FilterGenotypes","description":"Recalibrates qualities and filters genotypes","source":"@site/docs/modules/filter_genotypes.md","sourceDirName":"modules","slug":"/modules/fg","permalink":"/gatk-sv/docs/modules/fg","draft":false,"unlisted":false,"editUrl":"https://github.com/broadinstitute/gatk-sv/tree/master/website/docs/modules/filter_genotypes.md","tags":[],"version":"current","sidebarPosition":19,"frontMatter":{"title":"FilterGenotypes","description":"Recalibrates qualities and filters genotypes","sidebar_position":19,"slug":"fg"},"sidebar":"tutorialSidebar","previous":{"title":"SVConcordance","permalink":"/gatk-sv/docs/modules/svc"},"next":{"title":"AnnotateVcf","permalink":"/gatk-sv/docs/modules/av"}}');var t=l(4848),r=l(8453),s=l(1944);const o={title:"FilterGenotypes",description:"Recalibrates qualities and filters genotypes",sidebar_position:19,slug:"fg"},a=void 0,d={},c=[{value:"Model features",id:"model-features",level:3},{value:"Model availability",id:"model-availability",level:3},{value:"SL scores",id:"sl-scores",level:3},{value:"Modes",id:"modes",level:3},{value:"QC recommendations",id:"qc-recommendations",level:3},{value:"Inputs",id:"inputs",level:3},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>vcf</code>",id:"optional-vcf",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>output_prefix</code>",id:"optional-output_prefix",level:4},{value:"<code>ploidy_table</code>",id:"ploidy_table",level:4},{value:"<code>gq_recalibrator_model_file</code>",id:"gq_recalibrator_model_file",level:4},{value:"<code>recalibrate_gq_args</code>",id:"recalibrate_gq_args",level:4},{value:"<code>genome_tracks</code>",id:"genome_tracks",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>no_call_rate_cutoff</code>",id:"optional-no_call_rate_cutoff",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>fmax_beta</code>",id:"optional-fmax_beta",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>truth_json</code>",id:"optional-truth_json",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>sl_filter_args</code>",id:"optional-sl_filter_args",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>run_qc</code>",id:"optional-run_qc",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>optimize_vcf_records_per_shard</code>",id:"optional-optimize_vcf_records_per_shard",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>filter_vcf_records_per_shard</code>",id:"optional-filter_vcf_records_per_shard",level:4},{value:"Outputs",id:"outputs",level:3},{value:"<code>filtered_vcf</code>",id:"filtered_vcf",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>main_vcf_qc_tarball</code>",id:"optional-main_vcf_qc_tarball",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>vcf_optimization_table</code>",id:"optional-vcf_optimization_table",level:4},{value:"<HighlightOptionalArg>Optional</HighlightOptionalArg> <code>sl_cutoff_qc_tarball</code>",id:"optional-sl_cutoff_qc_tarball",level:4},{value:"<code>unfiltered_recalibrated_vcf</code>",id:"unfiltered_recalibrated_vcf",level:4}];function h(e){const i={a:"a",code:"code",em:"em",h3:"h3",h4:"h4",li:"li",mermaid:"mermaid",ol:"ol",p:"p",pre:"pre",ul:"ul",...(0,r.R)(),...e.components};return(0,t.jsxs)(t.Fragment,{children:[(0,t.jsx)(i.p,{children:(0,t.jsx)(i.a,{href:"https://github.com/broadinstitute/gatk-sv/blob/main/wdl/FilterGenotypes.wdl",children:"WDL source code"})}),"\n",(0,t.jsxs)(i.p,{children:["Performs genotype quality recalibration using a machine learning model based on ",(0,t.jsx)(i.a,{href:"https://github.com/dmlc/xgboost",children:"xgboost"}),"\nand filters genotypes. The output VCF contains the following updated fields:"]}),"\n",(0,t.jsxs)(i.ul,{children:["\n",(0,t.jsxs)(i.li,{children:[(0,t.jsx)(i.code,{children:"SL"})," : Scaled logit scores (see ",(0,t.jsx)(i.a,{href:"#sl-scores",children:"here"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:[(0,t.jsx)(i.code,{children:"GQ"})," : Updated genotype quality rescaled using ",(0,t.jsx)(i.code,{children:"SL"})]}),"\n",(0,t.jsxs)(i.li,{children:[(0,t.jsx)(i.code,{children:"OGQ"})," : Original ",(0,t.jsx)(i.code,{children:"GQ"})," score before recalibration"]}),"\n",(0,t.jsxs)(i.li,{children:[(0,t.jsx)(i.code,{children:"HIGH_NCR"})," : Filter status assigned to variants exceeding a ",(0,t.jsx)(i.a,{href:"#optional-no_call_rate_cutoff",children:"threshold proportion"}),"\nof no-call genotypes. This will also be applied to variants with genotypes that have already been filtered in the input VCF."]}),"\n"]}),"\n",(0,t.jsx)(i.p,{children:"The following diagram illustrates the recommended invocation order:"}),"\n",(0,t.jsx)(i.mermaid,{value:"\nstateDiagram\n  direction LR\n  \n  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d\n  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white\n  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d\n\n  svc: SVConcordance\n  fg: FilterGenotypes\n  avcf: AnnotateVcf\n  svc --\x3e fg\n  fg --\x3e avcf\n  \n  class fg thisModule\n  class svc inModules\n  class avcf outModules"}),"\n",(0,t.jsx)(i.h3,{id:"model-features",children:"Model features"}),"\n",(0,t.jsx)(i.p,{children:"The model uses the following features:"}),"\n",(0,t.jsxs)(i.ul,{children:["\n",(0,t.jsxs)(i.li,{children:["Genotype properties:","\n",(0,t.jsxs)(i.ul,{children:["\n",(0,t.jsx)(i.li,{children:"Non-reference and no-call allele counts"}),"\n",(0,t.jsxs)(i.li,{children:["Genotype quality (",(0,t.jsx)(i.code,{children:"GQ"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Supporting evidence types (",(0,t.jsx)(i.code,{children:"EV"}),") and respective genotype qualities (",(0,t.jsx)(i.code,{children:"PE_GQ"}),", ",(0,t.jsx)(i.code,{children:"SR_GQ"}),", ",(0,t.jsx)(i.code,{children:"RD_GQ"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Raw call concordance (",(0,t.jsx)(i.code,{children:"CONC_ST"}),")"]}),"\n"]}),"\n"]}),"\n",(0,t.jsxs)(i.li,{children:["Variant properties:","\n",(0,t.jsxs)(i.ul,{children:["\n",(0,t.jsxs)(i.li,{children:["Variant type (",(0,t.jsx)(i.code,{children:"SVTYPE"}),") and size (",(0,t.jsx)(i.code,{children:"SVLEN"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Calling algorithms (",(0,t.jsx)(i.code,{children:"ALGORITHMS"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Supporting evidence types (",(0,t.jsx)(i.code,{children:"EVIDENCE"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Two-sided SR support flag (",(0,t.jsx)(i.code,{children:"BOTHSIDES_SUPPORT"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Evidence overdispersion flag (",(0,t.jsx)(i.code,{children:"PESR_GT_OVERDISPERSION"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["SR noise flag (",(0,t.jsx)(i.code,{children:"HIGH_SR_BACKGROUND"}),")"]}),"\n",(0,t.jsxs)(i.li,{children:["Raw call concordance (",(0,t.jsx)(i.code,{children:"STATUS"}),", ",(0,t.jsx)(i.code,{children:"NON_REF_GENOTYPE_CONCORDANCE"}),", ",(0,t.jsx)(i.code,{children:"VAR_PPV"}),", ",(0,t.jsx)(i.code,{children:"VAR_SENSITIVITY"}),", ",(0,t.jsx)(i.code,{children:"TRUTH_AF"}),")"]}),"\n"]}),"\n"]}),"\n",(0,t.jsxs)(i.li,{children:["Reference context with respect to UCSC Genome Browser tracks:","\n",(0,t.jsxs)(i.ul,{children:["\n",(0,t.jsx)(i.li,{children:"RepeatMasker"}),"\n",(0,t.jsx)(i.li,{children:"Segmental duplications"}),"\n",(0,t.jsx)(i.li,{children:"Simple repeats"}),"\n",(0,t.jsx)(i.li,{children:"K-mer mappability (umap_s100 and umap_s24)"}),"\n"]}),"\n"]}),"\n"]}),"\n",(0,t.jsx)(i.h3,{id:"model-availability",children:"Model availability"}),"\n",(0,t.jsx)(i.p,{children:"For ease of use, we provide a model pre-trained on high-quality data with truth data derived from long-read calls:"}),"\n",(0,t.jsx)(i.pre,{children:(0,t.jsx)(i.code,{children:"gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/gatk-sv-recalibrator.aou_phase_1.v1.model\n"})}),"\n",(0,t.jsxs)(i.p,{children:['See the SV "Genotype Filter" section on page 34 of the ',(0,t.jsx)(i.a,{href:"https://support.researchallofus.org/hc/en-us/articles/4617899955092-All-of-Us-Genomic-Quality-Report-ARCHIVED-C2022Q4R9-CDR-v7",children:"All of Us Genomic Quality Report C2022Q4R9 CDR v7"})," for further details on model training. The generation and release of this model was made possible by the All of Us program (see ",(0,t.jsx)(i.a,{href:"/docs/acknowledgements",children:"here"}),")."]}),"\n",(0,t.jsx)(i.h3,{id:"sl-scores",children:"SL scores"}),"\n",(0,t.jsxs)(i.p,{children:['All valid genotypes are annotated with a "scaled logit" (',(0,t.jsx)(i.code,{children:"SL"}),") score, which is rescaled to non-negative adjusted ",(0,t.jsx)(i.code,{children:"GQ"})," values on [1, 99]. Note that the rescaled ",(0,t.jsx)(i.code,{children:"GQ"})," values should ",(0,t.jsx)(i.em,{children:"not"})," be interpreted as probabilities. Original genotype qualities are retained in the ",(0,t.jsx)(i.code,{children:"OGQ"})," field."]}),"\n",(0,t.jsxs)(i.p,{children:["A more positive ",(0,t.jsx)(i.code,{children:"SL"})," score indicates higher probability that the given genotype is not homozygous for the reference allele. Genotypes are therefore filtered using ",(0,t.jsx)(i.code,{children:"SL"})," thresholds that depend on SV type and size. This workflow also generates QC plots using the ",(0,t.jsx)(i.a,{href:"https://github.com/broadinstitute/gatk-sv/blob/main/wdl/MainVcfQc.wdl",children:"MainVcfQc"})," workflow to review call set quality (see below for recommended practices)."]}),"\n",(0,t.jsx)(i.h3,{id:"modes",children:"Modes"}),"\n",(0,t.jsx)(i.p,{children:"This workflow can be run in one of two modes:"}),"\n",(0,t.jsxs)(i.ol,{children:["\n",(0,t.jsxs)(i.li,{children:["\n",(0,t.jsxs)(i.p,{children:["(Recommended) The user explicitly provides a set of ",(0,t.jsx)(i.code,{children:"SL"})," cutoffs through the ",(0,t.jsx)(i.code,{children:"sl_filter_args"})," parameter, e.g."]}),"\n",(0,t.jsx)(i.pre,{children:(0,t.jsx)(i.code,{children:'"--small-del-threshold 93 --medium-del-threshold 150 --small-dup-threshold -51 --medium-dup-threshold -4 --ins-threshold -13 --inv-threshold -19"\n'})}),"\n",(0,t.jsxs)(i.p,{children:["Genotypes with ",(0,t.jsx)(i.code,{children:"SL"})," scores less than the cutoffs are set to no-call (",(0,t.jsx)(i.code,{children:"./."}),"). The above values were taken directly from Appendix N of the ",(0,t.jsx)(i.a,{href:"https://support.researchallofus.org/hc/en-us/articles/4617899955092-All-of-Us-Genomic-Quality-Report-ARCHIVED-C2022Q4R9-CDR-v7",children:"All of Us Genomic Quality Report C2022Q4R9 CDR v7 "}),". Users should adjust the thresholds depending on data quality and desired accuracy. Please see the arguments in ",(0,t.jsx)(i.a,{href:"https://github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/apply_sl_filter.py",children:"this script"})," for all available options."]}),"\n"]}),"\n",(0,t.jsxs)(i.li,{children:["\n",(0,t.jsxs)(i.p,{children:["(Advanced) The user provides truth labels for a subset of non-reference calls, and ",(0,t.jsx)(i.code,{children:"SL"})," cutoffs are automatically optimized. These truth labels should be provided as a json file in the following format:"]}),"\n",(0,t.jsx)(i.pre,{children:(0,t.jsx)(i.code,{className:"language-json",children:' {\n   "sample_1": \n   {\n     "good_variant_ids": ["variant_1", "variant_3"], \n     "bad_variant_ids": ["variant_5", "variant_10"]\n   },\n   "sample_2": \n   {\n     "good_variant_ids": ["variant_2", "variant_13"], \n     "bad_variant_ids": ["variant_8", "variant_11"]\n   }\n }\n'})}),"\n",(0,t.jsxs)(i.p,{children:['where "good_variant_ids" and "bad_variant_ids" are lists of variant IDs corresponding to non-reference (i.e. het or hom-var) sample genotypes that are true positives and false positives, respectively. ',(0,t.jsx)(i.code,{children:"SL"})," cutoffs are optimized by maximizing the ",(0,t.jsx)(i.a,{href:"https://en.wikipedia.org/wiki/F-score",children:"F-score"}),' with "beta" parameter ',(0,t.jsx)(i.code,{children:"fmax_beta"}),", which modulates the weight given to precision over recall (lower values give higher precision)."]}),"\n"]}),"\n"]}),"\n",(0,t.jsxs)(i.p,{children:['In both modes, the workflow additionally filters variants based on the "no-call rate", the proportion of genotypes that were filtered in a given variant. Variants exceeding the ',(0,t.jsx)(i.code,{children:"no_call_rate_cutoff"})," are assigned a ",(0,t.jsx)(i.code,{children:"HIGH_NCR"})," filter status."]}),"\n",(0,t.jsx)(i.h3,{id:"qc-recommendations",children:"QC recommendations"}),"\n",(0,t.jsxs)(i.p,{children:["We strongly recommend performing call set QC after this module. By default, QC plotting is enabled with the ",(0,t.jsx)(i.a,{href:"#optional-run_qc",children:"run_qc"}),"\nargument. Users should carefully inspect the main plots from the ",(0,t.jsx)(i.a,{href:"#optional-main_vcf_qc_tarball",children:"main_vcf_qc_tarball"}),".\nPlease see the ",(0,t.jsx)(i.a,{href:"./mvqc",children:"MainVcfQc"})," module documentation for more information on interpreting these plots and recommended\nQC criteria."]}),"\n",(0,t.jsx)(i.h3,{id:"inputs",children:"Inputs"}),"\n",(0,t.jsxs)(i.h4,{id:"optional-vcf",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"vcf"})]}),"\n",(0,t.jsxs)(i.p,{children:["Input VCF generated from ",(0,t.jsx)(i.a,{href:"./svc#concordance_vcf",children:"SVConcordance"}),"."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-output_prefix",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"output_prefix"})]}),"\n",(0,t.jsx)(i.p,{children:"Default: use input VCF filename. Prefix for the output VCF, such as the cohort name. May be alphanumeric with underscores."}),"\n",(0,t.jsx)(i.h4,{id:"ploidy_table",children:(0,t.jsx)(i.code,{children:"ploidy_table"})}),"\n",(0,t.jsxs)(i.p,{children:["Table of sample ploidies generated in ",(0,t.jsx)(i.a,{href:"./jrc#ploidy_table",children:"JoinRawCalls"}),"."]}),"\n",(0,t.jsx)(i.h4,{id:"gq_recalibrator_model_file",children:(0,t.jsx)(i.code,{children:"gq_recalibrator_model_file"})}),"\n",(0,t.jsxs)(i.p,{children:["GQ-Recalibrator model. A public model is listed as ",(0,t.jsx)(i.code,{children:"aou_recalibrate_gq_model_file"})," ",(0,t.jsx)(i.a,{href:"https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json",children:"here"}),"."]}),"\n",(0,t.jsx)(i.h4,{id:"recalibrate_gq_args",children:(0,t.jsx)(i.code,{children:"recalibrate_gq_args"})}),"\n",(0,t.jsxs)(i.p,{children:["Arguments to pass to the ",(0,t.jsx)(i.code,{children:"GQ"})," recalibration tool. Users should leave this with the default configuration in Terra."]}),"\n",(0,t.jsx)(i.h4,{id:"genome_tracks",children:(0,t.jsx)(i.code,{children:"genome_tracks"})}),"\n",(0,t.jsx)(i.p,{children:"Genome tracks for sequence context annotation. Users should leave this with the default configuration in Terra."}),"\n",(0,t.jsxs)(i.h4,{id:"optional-no_call_rate_cutoff",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"no_call_rate_cutoff"})]}),"\n",(0,t.jsxs)(i.p,{children:["Default: ",(0,t.jsx)(i.code,{children:"0.05"}),". Threshold fraction of samples that must have no-call genotypes in order to filter a variant. Set to 1 to disable."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-fmax_beta",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"fmax_beta"})]}),"\n",(0,t.jsxs)(i.p,{children:["Default: ",(0,t.jsx)(i.code,{children:"0.4"}),". If providing a truth set, defines the beta parameter for F-score optimization."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-truth_json",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"truth_json"})]}),"\n",(0,t.jsxs)(i.p,{children:["Truth labels for input variants. If provided, the workflow will attempt to optimize filtering cutoffs automatically\nusing the F-score. If provided, ",(0,t.jsx)(i.a,{href:"#optional-sl_filter_args",children:"sl_filter_args"})," is ignored."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-sl_filter_args",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"sl_filter_args"})]}),"\n",(0,t.jsxs)(i.p,{children:["Arguments for the ",(0,t.jsx)(i.a,{href:"https://github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/scripts/apply_sl_filter.py",children:"SL filtering script"}),".\nThis should be used to set ",(0,t.jsx)(i.code,{children:"SL"})," cutoffs for filtering (refer to description above). Overridden by ",(0,t.jsx)(i.a,{href:"#optional-truth_json",children:"truth_json"}),"."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-run_qc",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"run_qc"})]}),"\n",(0,t.jsxs)(i.p,{children:["Default: ",(0,t.jsx)(i.code,{children:"true"}),". Enable running ",(0,t.jsx)(i.a,{href:"./mvqc",children:"MainVcfQc"})," automatically. By default, filtered variants will be excluded from\nthe plots."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-optimize_vcf_records_per_shard",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"optimize_vcf_records_per_shard"})]}),"\n",(0,t.jsxs)(i.p,{children:["Default: ",(0,t.jsx)(i.code,{children:"50000"}),". Shard size for scattered cutoff optimization tasks. Decrease this if those steps are running slowly."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-filter_vcf_records_per_shard",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"filter_vcf_records_per_shard"})]}),"\n",(0,t.jsxs)(i.p,{children:["Default: ",(0,t.jsx)(i.code,{children:"20000"}),". Shard size for scattered ",(0,t.jsx)(i.code,{children:"GQ"})," recalibration tasks. Decrease this if those steps are running slowly."]}),"\n",(0,t.jsx)(i.h3,{id:"outputs",children:"Outputs"}),"\n",(0,t.jsx)(i.h4,{id:"filtered_vcf",children:(0,t.jsx)(i.code,{children:"filtered_vcf"})}),"\n",(0,t.jsx)(i.p,{children:"Filtered VCF."}),"\n",(0,t.jsxs)(i.h4,{id:"optional-main_vcf_qc_tarball",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"main_vcf_qc_tarball"})]}),"\n",(0,t.jsxs)(i.p,{children:["QC plots generated with ",(0,t.jsx)(i.a,{href:"./mvqc",children:"MainVcfQc"}),". Only generated if using ",(0,t.jsx)(i.a,{href:"#optional-run_qc",children:"run_qc"}),"."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-vcf_optimization_table",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"vcf_optimization_table"})]}),"\n",(0,t.jsxs)(i.p,{children:["Table of cutoff optimization metrics. Only generated if ",(0,t.jsx)(i.a,{href:"#optional-truth_json",children:"truth_json"})," is provided."]}),"\n",(0,t.jsxs)(i.h4,{id:"optional-sl_cutoff_qc_tarball",children:[(0,t.jsx)(s.$,{children:"Optional"})," ",(0,t.jsx)(i.code,{children:"sl_cutoff_qc_tarball"})]}),"\n",(0,t.jsxs)(i.p,{children:["Cutoff optimization and QC plots. Only generated if ",(0,t.jsx)(i.a,{href:"#optional-truth_json",children:"truth_json"})," is provided."]}),"\n",(0,t.jsx)(i.h4,{id:"unfiltered_recalibrated_vcf",children:(0,t.jsx)(i.code,{children:"unfiltered_recalibrated_vcf"})}),"\n",(0,t.jsxs)(i.p,{children:["Supplemental output of the VCF after assigning ",(0,t.jsx)(i.code,{children:"SL"})," genotype scores but before applying filtering."]})]})}function p(e={}){const{wrapper:i}={...(0,r.R)(),...e.components};return i?(0,t.jsx)(i,{...e,children:(0,t.jsx)(h,{...e})}):h(e)}},1944:(e,i,l)=>{l.d(i,{$:()=>t});var n=l(4848);const t=e=>{let{children:i}=e;return(0,n.jsx)("span",{style:{backgroundColor:"var(--highlight-optional-arg-background-color)",borderRadius:"2px",color:"var(--highlight-optional-arg-text-color)",padding:"0.2rem"},children:i})}},8453:(e,i,l)=>{l.d(i,{R:()=>s,x:()=>o});var n=l(6540);const t={},r=n.createContext(t);function s(e){const i=n.useContext(r);return n.useMemo((function(){return"function"==typeof e?e(i):{...i,...e}}),[i,e])}function o(e){let i;return i=e.disableParentContext?"function"==typeof e.components?e.components(t):e.components||t:s(e.components),n.createElement(r.Provider,{value:i},e.children)}}}]);