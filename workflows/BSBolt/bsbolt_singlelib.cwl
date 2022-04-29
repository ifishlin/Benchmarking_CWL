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
  DB:
    type: Directory
  threads:
    type: int
  output_name:
    type: string
  pbat:
    type: boolean
    default: False

steps:
  bsbolt_align:
     run: "./tools/bsbolt_align.cwl"
     scatter: [r1, r2]
     scatterMethod: 'dotproduct'
     in:
        r1:
          source: r1
        r2:
          source: r2
        DB:
          source: DB
        threads:
          source: threads
        output_name:
          source: output_name
        pbat:
          source: pbat
     out:
        - bam

  samtools_fixmate:
    run: "../../tools/samtools_fixmate.cwl"
    scatter: [bam]
    scatterMethod: 'dotproduct'
    in:
      bam: 
        source: bsbolt_align/bam
    out:
      - bam_fixmate

  merge_and_sort:
    run: "../../tools/samtools_merge_and_sort.cwl"
    in:
      bams:
        source: samtools_fixmate/bam_fixmate
      name_sort:
        valueFrom: $(false)
      threads: threads
    out:
      - bam_merged

  samtools_markdup:
    run: "../../tools/samtools_markdup.cwl"
    in:
      bam_sorted:
        source: merge_and_sort/bam_merged
    out:
      - bam_markdup

  samtools_index:
     run: "../../tools/samtools_index.cwl"
     in:
       bam_sorted:
         source: samtools_markdup/bam_markdup
     out:
       - bam_sorted_indexed

  bsbolt_call:
     run: "./tools/bsbolt_callMeth.cwl"
     in:
       DB:
         source: DB
       threads:
         source: threads
       output_name:
         source: output_name
       bam:
         source: samtools_index/bam_sorted_indexed
     out:
       - bggz

outputs: 
#  bam_fixmate:
#    type: File[]
#    outputSource: samtools_fixmate/bam_fixmate
   bam_merged:
     type: File
     outputSource: merge_and_sort/bam_merged
   bam_sorted_indexed:
     type: File
     outputSource: samtools_index/bam_sorted_indexed
   bggz:
     type: File
     outputSource: bsbolt_call/bggz
