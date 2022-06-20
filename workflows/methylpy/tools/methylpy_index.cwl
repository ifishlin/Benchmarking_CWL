#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["methylpy", "build-reference"]
arguments:
#  - valueFrom: "True"
#    prefix: "--bowtie2"
#    position: 3
  - valueFrom: $(inupts.output_name)
    prefix: "--output-prefix"
    position: 2

requirements:
  DockerRequirement:
    dockerPull: ifishlin324/methylpy
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
      prefix: "--input-files"
  - id: output_name
    type: string

outputs: 
  idx:
    type: File
    secondaryFiles:
      - ^_f.1.bt2
      - ^_f.2.bt2
      - ^_f.3.bt2
      - ^_f.4.bt2
      - ^_r.1.bt2
      - ^_r.2.bt2
      - ^_r.3.bt2
      - ^_r.4.bt2
      - ^_r.rev.1.bt2
      - ^_r.rev.2.bt2
      - ^_f.rev.1.bt2
      - ^_f.rev.2.bt2
    outputBinding:
      glob: $(inupts.output_name)
