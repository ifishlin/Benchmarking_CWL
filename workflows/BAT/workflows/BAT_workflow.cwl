#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
inputs:
  - id : reference
    type: File
  - id: query
    type: 'File[]'
  - id: mate_pair
    type: 'File[]'
  - id: prefix_db
    type: string
  - id: prefix_location
    type: Directory
  - id: path_outfiles
    type: string

outputs:
  bam:
    type: File
    outputSource: mapping/bam
  vcf_file:
    type: File
    outputSource: calling/vcf

steps:
  mapping:
    run: "../tools/BAT_mapping.cwl"
    in:
      reference: reference
      query: query
      mate_pair: mate_pair
      prefix_db: prefix_db
      prefix_location: prefix_location
      path_outfiles: path_outfiles
    out: [bam]

  calling:
    run: "../tools/BAT_calling.cwl"
    in: 
      reference: reference
      query: mapping/bam
    out: [vcf]
