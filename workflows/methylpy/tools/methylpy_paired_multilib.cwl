#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
 StepInputExpressionRequirement: {}
 InlineJavascriptRequirement: {}

inputs:
  read1:
    doc: first reads belonging to the same library
    type:
      type: array
      items: File
  read2:
    doc: first reads belonging to the same library
    type:
      type: array
      items: File
  sample:
    type: string
  ref_fasta:
    type: File
    secondaryFiles:
      - _f.1.bt2
      - _f.2.bt2
      - _f.3.bt2
      - _f.4.bt2
      - _f.rev.1.bt2
      - _f.rev.2.bt2
      - _r.1.bt2
      - _r.2.bt2
      - _r.3.bt2
      - _r.4.bt2
      - _r.rev.1.bt2
      - _r.rev.2.bt2
  pbat:
     type: boolean

steps: 
  methylpy_paired:
    run: "./methylpy_paired.cwl"
    in:
       read1:
         source: read1
       read2:
         source: read2
       sample:
         source: sample
       ref_fasta:
         source: ref_fasta
       pbat:
         valueFrom: $(false)
    out:
       - tsvgz

outputs: 
   tsvgz:
     type: File
     outputSource: methylpy_paired/tsvgz
