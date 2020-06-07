#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: Samtools process for unmapped reads aligned by Bowtie2 
doc: |
  This CWL converts aligned SAM files to BAM using samtools-view. Samtools-reheader is used 
  as Bowtie2 returns different header than Hisat2. Finally, satools-sort is used to acquire the final 
  unmapped - remapped BAM.
  version: 0.01
  
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  cores:
    type: int?
  input_sam:
    type: File[]
  header:
    type: File[]
  output_as_bam:
    type: boolean?
    default: true
  input_as_sam:
    type: boolean?
    default: true
  include_header:
    type: boolean?
    default: false

outputs:
  remap_bam:
    type: File[]
    outputSource: samtools_sort/output_sorted

steps:
  samtools_view_bam:
    run: ../tools/samtools-view.cwl
    scatter: input_sam
    in:
      output_as_bam: output_as_bam
      input_as_sam: input_as_sam
      include_header: include_header
      input_sam: input_sam
      bam_output:
        valueFrom: $(inputs.input_sam.nameroot).bam
    out: [output_bam]

  samtools_reheader:
    run: ../tools/samtools-reheader.cwl
    scatter:
      - input_bam
      - input_header
    scatterMethod: dotproduct
    in:
      input_bam: samtools_view_bam/output_bam
      input_header: header
      bam_output:
        valueFrom: $(inputs.input_bam.nameroot).uns
    out: [output_bam]

  samtools_sort:
    run: ../tools/samtools-sort.cwl
    scatter: input_bam
    in:
      input_bam: samtools_reheader/output_bam
      sorted_bam:
        valueFrom: $(inputs.input_bam.nameroot).bam
    out: [output_sorted]

# Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl
