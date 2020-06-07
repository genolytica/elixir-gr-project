cwlVersion: v1.0
class: CommandLineTool

label: samtools index
doc: |
  This CWL executes samtools index for one BAM file.
  version: 0.01

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.input_bam) ]
  DockerRequirement:
    dockerPull: "biocontainers/samtools:v1.7.0_cv2"

baseCommand: ["samtools","index"]

inputs:
  cores:
    type: int?
    doc: "number of cores"
    inputBinding:
      position: 1
      prefix: -@
  input_bam:
    type: File?
    doc: "removed duplicates BAM"
    inputBinding:
      position: 2
      valueFrom: $(self.basename)
  
outputs:
  output_index:
    type: File
    doc: "BAM index"
    secondaryFiles: 
      - ".bai"
    outputBinding:
      glob: $(inputs.input_bam.basename)

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
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl