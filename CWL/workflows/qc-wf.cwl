#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: FastQC parallel and MultiQC
doc: |
  This CWL scatters FastQC quality control over multiple processors, one for
  each file and creates the MultiQC report
  version: 0.01

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  threads:
    type: int?
  outdir:
    type: string
  fastq:
    type: File[]

steps:
  fastqc:
    doc: "FastQC - Quality Control for multiple fastq"
    run: fastqc-wf.cwl
    in:
      threads: threads
      fastq: fastq
      outdir: outdir
    out:
      - fastqc_zip 
      - fastqc_html
  multiqc:
    doc: "MultiQC - Report"
    run: ../tools/multiqc.cwl
    in:
      qc_files_array_of_array:
        - fastqc/fastqc_zip
        - fastqc/fastqc_html
    out:
      - multiqc_zip 
      - multiqc_html
  
outputs:
  fastqc_zip:
    type: File[]
    outputSource: fastqc/fastqc_zip
  fastqc_html:
    type: File[]
    outputSource: fastqc/fastqc_html
  multiqc_zip:
    type: File
    outputSource: multiqc/multiqc_zip
  multiqc_html:
    type: File
    outputSource: multiqc/multiqc_html

## Metadata
#$namespaces:
#  s: https://schema.org/
#  edam: http://edamontology.org/

#$schemas:
# - https://schema.org/version/latest/schema.rdf
# - http://edamontology.org/EDAM_1.18.owl
