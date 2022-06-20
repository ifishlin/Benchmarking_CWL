#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [BAT_mapping_sam] 
arguments:
  - valueFrom: $(inputs.r1.nameroot)
    position: 5
    prefix: -o

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: ifishlin324/bat 
    dockerOutputDirectory: /opt
#  InitialWorkDirRequirement:
#    listing:
#      - $(inputs.prefix_location)

inputs: 
  reference:
    type: File
    doc: path/filename of reference genome fasta
    inputBinding:
      prefix: -g
      position: 1
  r1:
    type: File
    doc: path/filename of query sequences
    inputBinding:
      prefix: -q
      position: 2
  r2:
    type: File?
    doc: path/filename of mate pair sequences
    inputBinding:
      prefix: -p   
      position: 3  
  prefix_db:
    type: File 
    doc: path/prefix of database indices
    inputBinding:
      prefix: -i  
      position: 4
      valueFrom: $(inputs.prefix_db.path.split('.').slice(0, -1).join('.'))
  threads:
    type: int
    inputBinding:
      prefix: -t
      position: 5 
#  prefix_location:
#    type: Directory
#  r1.nameroot:
#    type: string  
#    doc: path/prefix of outfiles   
#    inputBinding:
#      prefix: -o 
#      position: 5

## OUTPUT PART      
outputs:           
  bam:
    type: File
    outputBinding:
      glob: $(inputs.r1.nameroot + ".bam")
  bam.bai:
    type: File
    outputBinding:
      glob: $(inputs.r1.nameroot + ".bam.bai")      
  excluded_bam:
    type: File
    outputBinding:
      glob: $(inputs.r1.nameroot + ".excluded.bam") 
  excluded_bam.bai:
    type: File
    outputBinding:
      glob: $(inputs.r1.nameroot + ".excluded.bam.bai")
  log:
    type: File
    outputBinding:
      glob: "*.log" 
#  unmapped.gz:
#    type: File
#    outputBinding:
#      glob: $(inputs.r1.nameroot + ".unmapped.gz")             
