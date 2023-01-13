cwlVersion: v1.0
class: CommandLineTool

label: samtools merge
doc: |
  This CWL executes samtools merge for an array of BAM files.
  version: 0.02

requirements:
  DockerRequirement:
    dockerPull: "biocontainers/samtools:v1.7.0_cv2"
    
baseCommand: ["samtools","merge"]

inputs:
  cores:
    type: int?
    doc: "number of cores"
    inputBinding:
      position: 1
      prefix: -@
  force_overwrite:
    type: boolean?
    doc: "force overwrite output file"
    inputBinding:
      position: 2
      prefix: -f 
  use_header:
    type: File?
    doc: "use header from file"
    inputBinding:
      position: 3
      prefix: -h
  name_sorted:
    type: boolean?
    doc: "input is name sorted"
    inputBinding:
      position: 4
      prefix: -n
  output_as_sam:
    type: boolean?
    doc: "output in SAM"
    inputBinding:
      position: 5
      prefix: -u
  unique_rg:
    type: boolean?
    doc: "emit unique RG"
    inputBinding:
      position: 6
      prefix: -c
  unique_pg:
    type: boolean?
    doc: "emit unique PG"
    inputBinding:
      position: 7
      prefix: -p
  bam_output:
    type: string
    doc: "output BAM file name"
    inputBinding:
      position: 8
  input_bam:
    type: File[]
    doc: "input BAM files"
    inputBinding:
      position: 9
  
outputs:
  output_merged:
    type: File
    doc: "merged BAM output"
    outputBinding:
      glob: $(inputs.bam_output)

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