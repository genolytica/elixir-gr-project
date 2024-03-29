#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: Complete Mapping and Quality Control pipeline for Single-end data
doc: |
  A workflow which:
  i) Runs Hisat2 on each fastq while generating a fastq file containing unmapped reads 
  ii) Runs Bowtie2 using the --very-sensitive-local option, on the unmapped reads
  iii) A subworkflow is employed to prepare BAM files, corresponding to the 
      unmapped-remaped reads, to be merged with the mapped reads form Hisat2
  iv) After the merged BAMs are created, they are sorted and assigned a user defined 
       name (samle identifier)
  v) Corresponding index file for each sample is generated.
  vi) At the same time it performs quality control over the FASTQ using fastqc and assembles the MultiQC report
  
requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  MultipleInputFeatureRequirement: {}
  
inputs:
  # Input FASTQs
  fastq:
    type: File[]

  #Cores
  cores_hisat2:
    type: int?
  cores_bowtie2:
    type: int? 
  cores_samtools:
    type: int?
  cores_qc:
    type: int?

  # HISAT2 alignment
  add_chr:
    type: boolean?
  hs_idx_basedir:
    type: Directory
  hs_idx_basename:
    type: string
  sam_output:
    type: string?
    default: output.sam
  unmapped_fastq:
    type: string?
    default: unmapped.fastq

  # BOWTIE2 alignment - unmapped
  local:
    type: string?
  bt_idx_basedir:
    type: Directory
  bt_idx_basename:
    type: string
  bt_output:
    type: string?
    default: bowtie.sam

  # Samtools-View arguments
  output_as_bam:
    type: boolean?
    default: true
  input_as_sam:
    type: boolean?
    default: true
  include_header:
    type: boolean?
    default: false
  bam_output:
    type: string
    default: random
  header_only:
    type: boolean?
    default: false

  # MSI - arguments
  force_overwrite:
    type: boolean?
  sorted_bam:
    type: string[] 
  
outputs:
  fastqc_zip:
    type: File[]
    outputSource: qc/fastqc_zip
  fastqc_html:
    type: File[]
    outputSource: qc/fastqc_html
  multiqc_zip:
    type: File
    outputSource: qc/multiqc_zip
  multiqc_html:
    type: File
    outputSource: qc/multiqc_html
  hisat2_sam:
    type: File[]
    outputSource: hisat2/output_sam
  hisat2_unmpd:
    type: File[]
    outputSource: hisat2/output_unmapped
  samtools_header:
    type: File[]
    outputSource: samtools_view_header/output_bam
  samtools_bam:
    type: File[]
    outputSource: samtools_view_bam/output_bam
  bowtie_sam:
    type: File[]
    outputSource: bowtie2/bowtie2_sam
  remap_bam:
    type: File[]
    outputSource: unmapped_wf/remap_bam
  final_bam:
    type: File[]
    outputSource: msi_wf/final_bam
  final_bai:
    type: File[]
    outputSource: msi_wf/final_bai

steps:
  qc:
    run: ../workflows/qc-wf.cwl
    in:
      threads: cores_qc
      fastq: fastq
    out: [fastqc_zip,  fastqc_html, multiqc_zip, multiqc_html]

  hisat2:
    run: ../tools/hisat2_single.cwl
    scatter: fastq
    in:
      cores: cores_hisat2
      hisat2_idx_basedir: hs_idx_basedir
      hisat2_idx_basename: hs_idx_basename
      fastq: fastq
      sam_output:
        valueFrom: $(inputs.fastq.nameroot).hisat2.sam
      unmapped_fastq:
        valueFrom: $(inputs.fastq.nameroot)_unmapped.fastq
    out: [output_sam, output_unmapped]

  bowtie2:
    run: ../tools/bowtie2_single.cwl
    scatter: fq
    in:
      cores: cores_bowtie2
      local_preset: local
      bowtie2_idx_basedir: bt_idx_basedir
      bowtie2_idx_basename: bt_idx_basename
      fq: hisat2/output_unmapped
      sam_output:
        valueFrom: $(inputs.fq.nameroot)_remap_tmp.sam
    out: [bowtie2_sam]

  samtools_view_header:
    run: ../tools/samtools-view.cwl
    scatter: input_sam
    in:
      cores: cores_samtools
      header_only: header_only
      input_sam: hisat2/output_sam
      bam_output:
        valueFrom: $(inputs.input_sam.nameroot).header.sam
    out: [output_bam]

  samtools_view_bam:
    run: ../tools/samtools-view.cwl
    scatter: input_sam
    in:
      cores: cores_samtools
      output_as_bam: output_as_bam
      input_as_sam: input_as_sam
      include_header: include_header
      input_sam: hisat2/output_sam
      bam_output:
        valueFrom: $(inputs.input_sam.nameroot).bam
    out: [output_bam]

  unmapped_wf:
    run: ../workflows/unmapped.cwl
    in:
      cores: cores_samtools
      input_sam: bowtie2/bowtie2_sam
      header: samtools_view_header/output_bam
      output_as_bam: output_as_bam
      input_as_sam: input_as_sam
      include_header: include_header
    out: [remap_bam]

  msi_wf:
    run: ../workflows/merge-sort-index.cwl
    in:
      cores: cores_samtools
      force_overwrite: force_overwrite
      input_bam:
        - samtools_view_bam/output_bam
        - unmapped_wf/remap_bam
      sorted_bam: sorted_bam
    out: [final_bam, final_bai]

# Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
 - http://edamontology.org/EDAM_1.18.owl
