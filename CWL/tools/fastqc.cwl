#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: FastQC quality control
doc: |
  This CWL executes the established FastQC quality control pipeline for one
  file. Basic parameters are supported.
  version: 0.01

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/fastqc:0.11.9--0"
  InitialWorkDirRequirement:
    listing: [ $(inputs.fastq) ]
    
baseCommand: "fastqc"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 2048

inputs:
  threads:
    type: int?
    doc: "number of threads"
    inputBinding:
      position: 1
      prefix: -t
  outdir:
    type: string?
    doc: "output directory"
    inputBinding:
      position: 2
      prefix: -o
      valueFrom: '$(runtime.outdir)'
  fastq:
    type: File
    doc: "fastq reads"
    inputBinding:
      position: 3
  
outputs:
  fastqc_zip:
    type: File
    doc: "fastqc zip results"
    outputBinding:
      glob: "*_fastqc.zip"
  fastqc_html:
    type: File
    doc: "fastqc html report"
    outputBinding:
      glob: "*_fastqc.html"

# Metadata
s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4199-0333
    s:email: mailto:pmoulos@hybridstat.com
    s:name: Panagiotis Moulos

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
 - http://edamontology.org/EDAM_1.18.owl
