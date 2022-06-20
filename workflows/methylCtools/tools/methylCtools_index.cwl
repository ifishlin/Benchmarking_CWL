#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
 StepInputExpressionRequirement: {}
 InlineJavascriptRequirement: {}

inputs:
  ref:
    type: File
  output_name:
    type: string

steps:
  methylCtools_fapos:
    run: "./methylCtools_fapos.cwl"
    in:
       ref:
         source: ref
       output_name:
         source: output_name
    out:
       - gz

  methylCtools_tabix:
    run: "./methylCtools_tabix.cwl"
    in:
      posgz:
        source: methylCtools_fapos/gz
    out:
      - posgztbi

  methylCtools_faconv:
    run: "./methylCtools_faconv.cwl"
    in:
      ref:
        source: ref
      output_name:
        source: output_name
    out:
      - convfq

  methylCtools_bwa:
    run: "./methylCtools_bwa.cwl"
    in:
      convfq:
        source: methylCtools_faconv/convfq
    out:
      - fq


outputs: 
  tabix:
   type: File
   outputSource: methylCtools_tabix/posgztbi
  bwafq:
   type: File
   outputSource: methylCtools_bwa/fq 
