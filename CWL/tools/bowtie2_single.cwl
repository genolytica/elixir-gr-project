#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: BOWTIE2 alignment for single-end data
doc: |
  This CWL executes BOWTIE2 alignment {for unmapped fastq} given 1 FASTQ file
  version: 0.01

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0"
    
baseCommand: ["bowtie2"]
arguments:
  - prefix: -S
    valueFrom: $(runtime.outdir)/$(inputs.sam_output)
  - prefix: -x
    valueFrom: $(inputs.bowtie2_idx_basedir.path)/$(inputs.bowtie2_idx_basename)

inputs:
  cores:
    type: int?
    doc: "Number of cores"
    inputBinding:
      prefix: -p
  local_preset:
    type: string?
    doc: "Local mode"
    inputBinding:
      prefix: --local
  dovetail:
    type: boolean?
    default: false
    doc: "Dovetail option"
    inputBinding:
      prefix: --dovetail
  bowtie2_idx_basedir:
    label: "Path to the directory the index for the reference genome are in"
    doc: "Path to the directory the index for the index files, such as .1.ht2 / etc exist."
    type: Directory
  bowtie2_idx_basename:
    label: "Basename of the hisat2 index files"
    doc: "Basename of the hisat2 index files, not including extensions like .1.ht2"
    type: string 
  fq:
    type: File
    inputBinding:
      prefix: -U
  sam_output:
    type: string?
    default: bowtie.sam
    inputBinding:

outputs:
  bowtie2_sam:
#    format: edam:format_2573
    type: File
    doc: "alignments in SAM format"
    outputBinding:
      glob: $(inputs.sam_output)

#Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-https.rdf
  - http://edamontology.org/EDAM_1.18.owl