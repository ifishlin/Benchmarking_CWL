#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["bwa", "index"]
arguments:
  - valueFrom: bwtsw
    position: 1
    prefix: -a


requirements:
  DockerRequirement:
    dockerPull: ifishlin324/methylctools
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.convfq)

inputs:
  - id: convfq
    type: File
    inputBinding:
      position: 2

outputs: 
  fq:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .aa
    outputBinding:
      glob: $(inputs.convfq.basename)
