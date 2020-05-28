#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: HISAT2 alignment for paired-end data
doc: |
  This CWL executes HISAT2 garph-based alignment given 2 FASTQ files from paired-end
  sequencing. Inputs fastq1 and fastq2 can be comma-separated lists of files containing
  mate 1s and 2s respectively.
  version: 0.01

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/hisat2:2.2.0--py36he1b5a44_1"
#   dockerPull: "quay.io/biocontainers/hisat2:2.2.0--py37he1b5a44_0"

baseCommand: ["hisat2"]
arguments:
  - prefix: --un-conc
    valueFrom: $(runtime.outdir)/$(inputs.unmapped_fastq)
  - prefix: -S
    valueFrom: $(runtime.outdir)/$(inputs.sam_output)
  - prefix: -x
    valueFrom: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename)

inputs:
  cores:
    type: int?
    doc: "number of cores"
    inputBinding:
      prefix: -p
  add_chr:
    type: boolean?
    doc: "add chr to reference"
    inputBinding:
      prefix: --add-chrname
  hisat2_idx_basedir:
    label: "Path to the directory the index for the reference genome are in"
    doc: "Path to the directory the index for the index files, such as .1.ht2 / etc exist."
    type: Directory
  hisat2_idx_basename:
    label: "Basename of the hisat2 index files"
    doc: "Basename of the hisat2 index files, not including extensions like .1.ht2"
    type: string
  fastq1:
    type: File
    doc: "fastq reads"
    inputBinding:
      prefix: "-1"
  fastq2:
    type: File
    doc: "fastq mates"
    inputBinding:
      prefix: "-2"
  sam_output:
    type: string?
    default: output.sam
  unmapped_fastq:
    type: string?
    default: unmapped.fastq
  
outputs:
  output_sam:
    type: File
    doc: "alignments in SAM format"
    outputBinding:
      glob: $(inputs.sam_output)
  output_unmapped1:
    type: File 
    doc: "unmapped reads 1 in FASTQ format"
    outputBinding:
      glob: "*1.fastq"
  output_unmapped2:
    type: File
    doc: "unmapped reads 2 in FASTQ format"
    outputBinding:
      glob: "*2.fastq" 

#stdout: $(inputs.sam_output)

# Metadata
#$namespaces:
#  s: https://schema.org/
#  edam: http://edamontology.org/

#$schemas:
#  - https://schema.org/docs/schema_org_rdfa.html
#  - http://edamontology.org/EDAM_1.18.owl
