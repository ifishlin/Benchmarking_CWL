#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["methylCtools", "fapos"]
arguments:
  - valueFrom: "-"
    position: 2
  - valueFrom: "|"
    position: 3
    shellQuote: false
  - valueFrom: bgzip
    position: 4
  - valueFrom: ">"
    position: 5
  - valueFrom: $(inputs.output_name).gz
    position: 6

requirements:
  DockerRequirement:
    dockerPull: ifishlin324/methylctools
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  - id: ref
    type: File
    inputBinding:
      position: 1
  - id: output_name
    type: string

outputs: 
  gz:
    type: File
    outputBinding:
      glob: $(inputs.output_name).gz
