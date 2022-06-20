#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
 StepInputExpressionRequirement: {}
 InlineJavascriptRequirement: {}

inputs:
  r1:
    type: File[]
  r2:
    type: File[]
  threads:
    type: int
  output_name:
    type: string
  ref_conv_fa:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  ref_pos:
    type: File
    secondaryFiles:
      - .tbi

steps:
  methylCtools_fqconv:
    run: "./tools/methylCtools_fqconv.cwl"
    scatter: [r1, r2]
    scatterMethod: 'dotproduct'
    in:
       r1:
         source: r1
       r2:
         source: r2
       output_name:
         source: output_name
    out:
       - convfq

  methylCtools_align:
    run: "./tools/methylCtools_align.cwl"
    scatter: [read_conv_fq]
    scatterMethod: 'dotproduct'
    in:
      ref_conv_fa:
        source: ref_conv_fa
      read_conv_fq:
        source: methylCtools_fqconv/convfq
      output_name:
        source: output_name
      threads:
        source: threads
    out:
      - bam

  methylCtools_bconv:
    run: "./tools/methylCtools_bconv.cwl"
    scatter: [bam]
    scatterMethod: 'dotproduct'
    in:
      bam:
        source: methylCtools_align/bam
      output_name:
        source: output_name
    out:
      - convbam

  merge_and_sort:
    run: "../tools/samtools_merge_and_sort.cwl"
    in:
      bams:
        source: methylCtools_bconv/convbam
      name_sort:
        valueFrom: $(false)
      threads: threads
    out:
       - bam_merged

  picard_markdup:
    run: "../tools/picard_markdup.cwl"
    in:
      bam_sorted:
        source:  merge_and_sort/bam_merged
    out:
      - bam_duprem
      - picard_markdup_log
      - picard_markdup_stat

  samtools_index:
     run: "../tools/samtools_index.cwl"
     in:
       bam_sorted:
         source: picard_markdup/bam_duprem
     out:
       - bam_sorted_indexed

  methylCtools_bcall:
     run: "./tools/methylCtools_bcall.cwl"
     in:
        bam:
          source: samtools_index/bam_sorted_indexed
        ref_pos:
          source: ref_pos
        output_name:
          source: output_name
     out:
        - callgz

  tabix:
    run: "../tools/tabix.cwl"
    in:
      vcfgz:
        source: methylCtools_bcall/callgz
    out:
      - vcfgztbi

outputs: 
   bam:
     type: File[]
     outputSource: methylCtools_align/bam
   convbam:
     type: File[]
     outputSource: methylCtools_bconv/convbam
   bam_sorted_indexed:
     type: File
     outputSource: samtools_index/bam_sorted_indexed
   vcfgztbi:
     type: File
     outputSource: tabix/vcfgztbi
