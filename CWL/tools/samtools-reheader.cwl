cwlVersion: v1.0
class: CommandLineTool

label: samtools reheader
doc: |
  This CWL executes samtools reheader operations to replace the header in given BAM (.bad)
  with the header in header.sam [samtools-view] and writes the new file (.uns)
  version: 0.01

requirements:
  DockerRequirement:
    dockerPull: "biocontainers/samtools:v1.7.0_cv4"
    
baseCommand: ["samtools","reheader"]

inputs:
  no_pg:
    type: boolean?
    default: false
    doc: "Add a @PG line to the header of the output file"
    inputBinding:
      prefix: -P
  in_place:
    type: boolean?
    default: false
    doc: "Perform the header edit in-place, if possible"
    inputBinding:
      prefix: -i
  input_header:
    type: File
    doc: "header input"
    inputBinding:
      position: 1
  input_bam:
    type: File
    doc: "BAM input"
    inputBinding:
      position: 2
  bam_output:
    type: string


outputs:
  output_bam:
    type: File
    doc: "Replaced header - BAM"
    outputBinding:
      glob: $(inputs.bam_output)

stdout: $(inputs.bam_output)

## Metadata
#$namespaces:
#  s: https://schema.org/
#  edam: http://edamontology.org/
#
#$schemas:
#  - https://schema.org/docs/schema_org_rdfa.html
#  - http://edamontology.org/EDAM_1.18.owl