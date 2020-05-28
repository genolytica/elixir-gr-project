cwlVersion: v1.0
class: CommandLineTool

label: samtools sort
doc: |
  This CWL executes samtools sort operations for one BAM file. Basic parameters
  are supported.
  version: 0.02

requirements:
  DockerRequirement:
    dockerPull: "biocontainers/samtools:v1.7.0_cv2"
    
baseCommand: ["samtools","sort"]

inputs:
  cores:
    type: int?
    doc: "number of cores"
    inputBinding:
#      position: 1
      prefix: -@
  name_sort:
    type: boolean?
    doc: "name sort"
    inputBinding:
#      position: 2
      prefix: -n
  input_bam:
    type: File?
    doc: "unsorted BAM"
    inputBinding:
#      position: 3
  sorted_bam:
    type: string
    doc: "sorted BAM"
    inputBinding:
      prefix: -o
  
outputs:
  output_sorted:
    type: File
    doc: "sorted BAM output"
    outputBinding:
      glob: $(inputs.sorted_bam)

stdout: $(inputs.sorted_bam)

## Metadata
#s:author:
#  - class: s:Person
#    s:identifier: https://orcid.org/0000-0002-4199-0333
#    s:email: mailto:pmoulos@hybridstat.com
#    s:name: Panagiotis Moulos

#$namespaces:
#  s: https://schema.org/
#  edam: http://edamontology.org/

#$schemas:
# - https://schema.org/docs/schema_org_rdfa.html
# - http://edamontology.org/EDAM_1.18.owl