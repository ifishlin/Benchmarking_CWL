#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["gsnap", "--gunzip",  "-O", "-A", "sam"] 
arguments:
  - valueFrom: "|"
    position: 6
    shellQuote: false
  - valueFrom: samtools
    position: 7
  - valueFrom: view
    position: 8
  - prefix: "-@"
    valueFrom: $(inputs.threads)
    position: 9
  - prefix: "-o"
    valueFrom: $(inputs.output_name).bam
    position: 10
  - valueFrom: "-"
    position: 11

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: ifishlin324/gsnap
#    dockerOutputDirectory: /data

stdout: stderr
stderr: $(inputs.output_name + ".gsnap.aln.log")

inputs: 
  - id: ref
    type: string
    inputBinding:
      prefix: -d
      position: 1
  - id: r1
    doc: read1.fa.gz
    type: File
    inputBinding:
      position: 4
  - id: r2
    doc: read2.fa.gz
    type: File
    inputBinding:
      position: 5
  - id: output_name
    type: string
  - id: threads
    type: int
    inputBinding:
      prefix: -t
      position: 2
  - id: index_dir
    type: Directory
    inputBinding:
      prefix: -D
      position: 3
  - id: pbat
    type: boolean
    inputBinding:
      position: 1
      valueFrom: ${if(inputs.pbat) return "--mode=cmet-nonstranded"; else return "--mode=cmet-stranded"}


outputs:
  bam:
    type: File
    outputBinding:
      glob: "*.bam"
  log:
    type: stderr

## OUTPUT PART      
#outputs: 
#  sam:
#    type: stdout
#  log:
#    type: stderr
#
#stdout: $(inputs.output_name).sam
#stderr: $(inputs.output_name).log

