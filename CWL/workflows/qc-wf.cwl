#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: FastQC parallel
doc: |
  This CWL scatters FastQC quality control over multiple processors, one for
  each file and MultiQC report
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
    run: ../workflows/fastqc-wf.cwl
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

# Metadata
#s:author:
#  - class: s:Person
#    s:identifier: https://orcid.org/0000-0002-4199-0333
#    s:email: mailto:pmoulos@hybridstat.com
#    s:name: Panagiotis Moulos
#
#$namespaces:
#  s: https://schema.org/
#  edam: http://edamontology.org/

#$schemas:
# - https://schema.org/docs/schema_org_rdfa.html
# - http://edamontology.org/EDAM_1.18.owl
