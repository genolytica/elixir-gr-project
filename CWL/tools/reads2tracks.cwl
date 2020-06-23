#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: metaseqR2 for Signal Tracks
doc: |

requirements:
  DockerRequirement:
    dockerPull: "daphnelettos/metaseqr2:1.3"

baseCommand: "Rscript"

inputs:
  script:
    type: File
    default:
      class: File
      path: ../runSignalTracks.R
    inputBinding:
      position: 1
  path:
    type: Directory
    doc: "Path where all the BED/BAM files are placed"
    default:
      class: Directory
      path: ../example_data
    inputBinding:
      prefix: --path
      position: 2
  targets:
#    type:
#      type: record
#      fields:
#        sampleid:
#          type: string[]
#        fileid:
#          type: string[]
#        condition:
#          type: string[]
#        reads:
#          type: string[]?
#        stranded:
#          type: string[]?
#   type:
#      type: array
#     items: string[]
#    type: 
#      type: array
#      items:
#        type: array
#        items: string
#    type: Any
    type: File
    doc: "targets file"
    inputBinding:
      prefix: --targets
      position: 3
  organism:
    type: string
    doc: "organism annotation"
    inputBinding:
      prefix: --org
      position: 4
  urlbase:
    type: string?
    doc: "A valid URL which is prepended to the created bigWig files"
    inputBinding:
      prefix: --urlbase
      position: 5
  stranded:
    type: boolean?
    doc: "Separate + and - strands and create separate bigWig files"
    inputBinding:
      prefix: --stranded
      position: 6
  normto:
    type: int?
    doc: "sum of signal for normalization"
    inputBinding:
      prefix: --normto
      position: 7
  exportpath:
    type: Directory?
    doc: "Path to export tracks"
    inputBinding:
      prefix: --exportpath
      position: 8
  hubname:
    type: string?
    doc: "Name of the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_name
      position: 9
  hubslbl:
    type: string?
    doc: "Short label for the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_sl
      position: 10
  hubllbl:
    type: string?
    doc: "Long label for the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_ll
      position: 11
  hubmail:
    type: string?
    doc: "Email associated with the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_email
      position: 12
  overwrite:
    type: boolean?
    doc: "Overwrite tracks if they exist"
    inputBinding:
      prefix: --overwrite
      position: 13
  rc:
    type: int?
    doc: "Fraction of cores to use"
    inputBinding:
      prefix: --rc
      position: 14

outputs:
  tracks_output:
    type: File[]
    outputBinding:
      glob: "*.txt"
  bigwig:
    type: File[]
    doc: "BigWig Files"
    outputBinding:
      glob: "*.bigWig"
  bigwig2:
    type: Directory?
    outputBinding:
      glob: $(inputs.organism)

#Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf
  - http://edamontology.org/EDAM_1.18.owl