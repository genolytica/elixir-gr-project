cwlVersion: v1.0
class: CommandLineTool

label: samtools view
doc: |
  This CWL executes samtools view operations for one BAM file. Basic parameters
  are supported.
  version: 0.03

requirements:
  DockerRequirement:
    dockerPull: "biocontainers/samtools:v1.7.0_cv4"
    
baseCommand: ["samtools","view"]

inputs:
  cores:
    type: int?
    doc: "number of cores"
    inputBinding:
      prefix: -@
  output_as_bam:
    type: boolean?
    doc: "output in BAM"
    inputBinding:
      position: 1
      prefix: -b
  input_as_sam:
    type: boolean?
    doc: "input in SAM"
    inputBinding:
      position: 2
      prefix: -S
  filter_exclude:
    type: int?
    doc: "do not output alignments not passing this filter"
    inputBinding:
      position: 3
      prefix: -F
  filter_include:
    type: int?
    doc: "do output alignments passing this filter"
    inputBinding:
      position: 4
      prefix: -f
  mapping_quality:
    type: int?
    doc: "skip alignment with mapping quality below this value"
    inputBinding:
      position: 5
      prefix: -q
  header_only:
    type: boolean?
    default: false
    doc: "output only the BAM header"
    inputBinding:
      position: 6
      prefix: -H
  count_only:
    type: boolean?
    default: false
    doc: "output only the number of BAM records"
    inputBinding:
      position: 7
      prefix: -c
  include_header:
    type: boolean?
    default: false
    doc: "include header"
    inputBinding:
      position: 8
      prefix: -h 
  input_sam:
    type: File?
    doc: "SAM input"
    inputBinding:
      position: 9
  bam_output:
    type: string
  
outputs:
  output_bam:
    type: File
    doc: "unsorted BAM output"
    outputBinding:
      glob: $(inputs.bam_output)

stdout: $(inputs.bam_output)

## Metadata
#s:author:
#  - class: s:Person
#    s:identifier: https://orcid.org/0000-0002-4199-0333
#    s:email: mailto:pmoulos@hybridstat.com
#    s:name: Panagiotis Moulos
#
#$namespaces:
#  s: https://schema.org/
#  edam: http://edamontology.org/
#
#$schemas:
# - https://schema.org/docs/schema_org_rdfa.html
# - http://edamontology.org/EDAM_1.18.owl