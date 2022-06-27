cwlVersion: v1.0
class: Workflow

requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
 StepInputExpressionRequirement: {}
 InlineJavascriptRequirement: {}

inputs:
  read1:
    type: File
  read2:
    type: File
  ref:
    type: File
    secondaryFiles:
      - _strands
  output_name:
    type: string
  pbat:
    type: boolean
    default: $(false)

steps:
  fame:
    run: "./tools/fame.cwl"
    in:
       read1:
         source: read1
       read2:
         source: read2
       ref:
         source: ref
       output_name:
         source: output_name
       pbat:
         source: pbat
    out:
       - tsv
       - log

outputs: 
   tsv:
     type: File
     outputSource: fame/tsv
   log:
     type: File
     outputSource: fame/log     
