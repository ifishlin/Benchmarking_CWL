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
  ref:
    type: File
    secondaryFiles:
      - .bis.amb
      - .bis.ann
      - .bis.pac
      - .dau.bwt
      - .dau.sa
      - .par.bwt
      - .par.sa
      - .fai
  pbat:
    type: boolean
    default: $(false)
  threads:
    type: int
  output_name:
    type: string 

steps:
  biscuit_align:
    run: "./tools/biscuit_align.cwl"
    scatter: [r1, r2]
    scatterMethod: 'dotproduct'
    in:
       r1:
         source: r1
       r2:
         source: r2
       ref:
         source: ref
       pbat:
         source: pbat
       threads:
         source: threads
       output_name:
         source: output_name
    out:
       - bam

  samblaster_sort:
    run: "../../tools/samblaster_sort.cwl"
    scatter: [bam]
    scatterMethod: 'dotproduct'
    in:
      bam: 
        source: biscuit_align/bam
      output_name:
        source: output_name
      threads:
        source: threads
    out:
      - bam_duprem
      - log

  samtools_merge:
    run: "../../tools/samtools_merge.cwl"
    in:
      bams:
        source: samblaster_sort/bam_duprem
      output_name:
        source: output_name
    out:
       - bam_merged

  samtools_index:
     run: "../../tools/samtools_index.cwl"
     in:
       bam_sorted:
         source: samtools_merge/bam_merged
     out:
       - bam_sorted_indexed

  biscuit_pileup:
     run: "./tools/biscuit_pileup.cwl"
     in:
       output_name: 
         source: output_name
       ref:
         source: ref
       bam_sorted:
         source: samtools_index/bam_sorted_indexed
     out:
       - vcf

  bgzip:
    run: "../../tools/bgzip.cwl"
    in:
      vcf:
        source: biscuit_pileup/vcf
    out:
      - vcf.gz

  tabix:
    run: "../../tools/tabix.cwl"
    in:
      vcfgz: 
        source: bgzip/vcf.gz
    out:
      - vcfgztbi

  biscuit_vcf2bed:
     run: "./tools/biscuit_vcf2bed.cwl"
     in:
       vcfgz:
         source: tabix/vcfgztbi
     out:
       - bed

outputs:
    samblaster_log:
      type: File[]
      outputSource: samblaster_sort/log 
    bams:
      type: File
      outputSource: samtools_index/bam_sorted_indexed
    vcfgztbi:
      type: File
      outputSource: tabix/vcfgztbi
    bed:
      type: File
      outputSource: biscuit_vcf2bed/bed