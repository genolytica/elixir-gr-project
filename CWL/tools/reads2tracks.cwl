#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: metaseqR2 for Signal Tracks
doc: |

requirements:
  DockerRequirement:
    dockerPull: "daphnelettos/metaseqr2:1.3"
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.path)
        writable: true
      - entry: $(inputs.exportpath)
        writable: true

baseCommand: "Rscript"

inputs:
  script:
    type: File?
    default:
      class: File
      path: ../runSignalTracks.R
    inputBinding:
      position: 1
  path:
    type: Directory
    doc: "Path where all the BED/BAM files are placed"
    inputBinding:
      prefix: --path
      position: 2
  samples:
    type: string
    doc: "Sample IDs - Enter as comma-separated list (no space)."
    inputBinding:
      prefix: --samples
      position: 3
  files:
    type: string
    doc: "File names - Enter as comma-separated list (no space)."
    inputBinding:
      prefix: --files
      position: 4
  conditions:
    type: string
    doc: "Sample conditions - Enter as comma-separated list (no space)."
    inputBinding:
      prefix: --conditions
      position: 5
  paired:
    type: string?
    doc: "Paired or single-end reads. Enter as comma-separated list (no space)."
    inputBinding:
      prefix: --paired
      position: 6
  strandp:
    type: string?
    doc: "Strand library construction protocol. Enter as comma-separated list (no space)."
    inputBinding:
      prefix: --strandp
      position: 7
#  targets:
#    type: File?
#    doc: "File array where BED/BAM files are placed"
#    inputBinding:
#      prefix: --targets
#      position: 8
  organism:
    type: string
    doc: "organism annotation"
    inputBinding:
      prefix: --org
      position: 9
  urlbase:
    type: string?
    doc: "A valid URL which is prepended to the created bigWig files"
    inputBinding:
      prefix: --urlbase
      position: 10
  stranded:
    type: boolean?
    doc: "Separate + and - strands and create separate bigWig files"
    inputBinding:
      prefix: --stranded
      position: 11
  normto:
    type: int?
    doc: "sum of signal for normalization"
    inputBinding:
      prefix: --normto
      position: 12
  exportpath:
    type: Directory?
    doc: "Path to export tracks"
    default:
      class: Directory
      path: ../output/signal_tracks
    inputBinding:
      prefix: --exportpath
      position: 13
  hubname:
    type: string?
    doc: "Name of the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_name
      position: 14
  hubslbl:
    type: string?
    doc: "Short label for the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_sl
      position: 15
  hubllbl:
    type: string?
    doc: "Long label for the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_ll
      position: 16
  hubmail:
    type: string?
    doc: "Email associated with the track hub created in case of stranded tracks"
    inputBinding:
      prefix: --hubinfo_email
      position: 17
  overwrite:
    type: boolean?
    doc: "Overwrite tracks if they exist"
    inputBinding:
      prefix: --overwrite
      position: 18
  rc:
    type: float?
    doc: "Fraction of cores to use"
    inputBinding:
      prefix: --rc
      position: 19

outputs:
  results:
    type:
      type: array
      items: [File, Directory]
    outputBinding:
      glob: $(inputs.exportpath.basename)


#Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf
  - http://edamontology.org/EDAM_1.18.owl
