#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: FastQC parallel
doc: |
  This CWL scatters FastQC quality control over multiple processors, one for
  each file.
  version: 0.11.9

requirements:
  ScatterFeatureRequirement: {}
    
inputs:
  threads:
    type: int?
  outdir:
    type: string?
  fastq:
    type: File[]

steps:
  fastqc:
    run: ../tools/fastqc.cwl
    scatter: fastq
    in:
      threads: threads
      fastq: fastq
      outdir: outdir
    out:
      - fastqc_zip 
      - fastqc_html
  
outputs:
  fastqc_zip:
    type: File[]
    outputSource: fastqc/fastqc_zip
  fastqc_html:
    type: File[]
    outputSource: fastqc/fastqc_html

# Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl