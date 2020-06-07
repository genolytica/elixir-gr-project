#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: Samtools Merge - Sort - Index
doc: |
  This CWL performs samtools-merge in order to combine aligned reads from Hisat2 and
  Bowtie2 (unmapped reads) for an array of file pairs. Subsequently the merged BAMs are
  sorted and named according to sample (user input) and the index file for every BAM
  is generated (.bai).
  version: 0.01

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  
inputs:
  cores:
    type: int?
  input_bam:
    type: 
      type: array
      items:
        type: array
        items: File
  force_overwrite:
    type: boolean?
  sorted_bam:
    type: string[]

outputs:
  final_bam:
    type: File[]
    outputSource: samtools_sort/output_sorted
  final_bai:
    type: File[]
    outputSource: samtools_index/output_index

steps:
  samtools_merge:
    run: ../tools/samtools-merge.cwl
    scatter: [input_bam, bam_output]
    scatterMethod: dotproduct
    in:
      cores: cores
      force_overwrite: force_overwrite
      bam_output: sorted_bam
      input_bam: input_bam
    out: [output_merged]

  samtools_sort:
    run: ../tools/samtools-sort.cwl
    scatter: [input_bam, sorted_bam]
    scatterMethod: dotproduct
    in:
      cores: cores
      input_bam: samtools_merge/output_merged
      sorted_bam: sorted_bam
    out: [output_sorted]

  samtools_index:
    run: ../tools/samtools-index.cwl
    scatter: input_bam
    in:
      cores: cores
      input_bam: samtools_sort/output_sorted
    out: [output_index]

# Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl
