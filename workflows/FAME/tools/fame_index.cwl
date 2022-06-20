#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["FAME"]

requirements:
  DockerRequirement:
    dockerPull: ifishlin324/fame
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.ref)

inputs:
  - id: ref
    type: File
    inputBinding:
      position: 1
      prefix: "--genome"
  - id: output_name
    type: string
    inputBinding:
      position: 2
      prefix: "--store_index"

outputs:
  genome:
    type: File
    secondaryFiles:
      - _strands
    outputBinding:
      glob: $(inputs.output_name)
