#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["biscuit", "align"]
arguments:
  - valueFrom: "|"
    position: 3
    shellQuote: false
  - valueFrom: samtools
    position: 4
  - valueFrom: sort
    position: 5
  - prefix: "-@"
    valueFrom: $(inputs.threads)
    position: 6
  - prefix: "-o"
    valueFrom: $(inputs.output_name)
    position: 7
  - valueFrom: "-"
    position: 8

hints:
  DockerRequirement:
    dockerPull: mgibio/biscuit
inputs:
  - id: r1
    type: File
    inputBinding:
      position: 4
  - id: r2
    type: File
    inputBinding:
      position: 5
  - id: ref
    type: Directory
    inputBinding:
      position: 3
  - id: direction
    type: int
    inputBinding:
      prefix: -b
      position: 1
  - id: threads
    type: int
    inputBinding:
      prefix: -t
      position: 2
  - id: output_name
    type: string
outputs:
  - id: bam:
    type: File
stdout: output.txt
