#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
 StepInputExpressionRequirement: {}
 InlineJavascriptRequirement: {}

inputs:
  - id : reference
    type: File
    secondaryFiles:
      - .fai
  - id: query
    type: File[]
  - id: mate_pair
    type: File[]
  - id: prefix_db
    type: File
    secondaryFiles:
      - ^.ctidx
      - ^.gaidx
  - id: threads
    type: int
  - id: output_name
    type: string
  - id: header
    type: File

outputs:
  bam_sorted_indexed:
    type: File
    outputSource: samtools_index/bam_sorted_indexed
  vcfgztbi:
    type: File
    outputSource: calling/vcfgztbi

steps:
  mapping:
    run: "./tools/BAT_mapping.cwl"
    scatter: [r1, r2]
    scatterMethod: 'dotproduct'
    in:
      reference: reference
      r1: query
      r2: mate_pair
      prefix_db: prefix_db
      threads: threads
    out: [bam]

  merge_and_sort:
    run: "../../tools/samtools_merge_and_sort.cwl"
    in:
      bams:
        source: mapping/bam
      name_sort:
        valueFrom: $(false)
      threads: threads
    out:
       - bam_merged

  samtools_index:
     run: "../../tools/samtools_index.cwl"
     in:
       bam_sorted:
         source: merge_and_sort/bam_merged
     out:
       - bam_sorted_indexed

  calling:
     run: "./tools/BAT_calling_latest.cwl"
     in:
       ref:
         source: reference
       bam:
         source: samtools_index/bam_sorted_indexed
       threads:
         source: threads
       output_name:
         source: output_name
       header:
         source: header
     out: [vcfgztbi]
