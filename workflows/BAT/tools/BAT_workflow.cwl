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
  - id: query
    type: File[]
  - id: mate_pair
    type: File[]
  - id: prefix_db
    type: File
    secondaryFiles:
      - .ctidx
      - .gaidx
#  - id: prefix_location
#    type: Directory
#  - id: path_outfiles
#    type: string

outputs:
  bam:
    type: File[]
    outputSource: mapping/bam

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
    out: [bam]
