#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: MultiQC quality control
doc: |
  This CWL executes the established MultiQC pipeline across multiple
  fastqc files to generate a single quality control report.
  version: 0.01
 
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/multiqc:1.8--py_2"
  InitialWorkDirRequirement:
    # This step is necessary since the input files
    # must be loaded into the working directory as there
    # is no way to specify the input file directly on the
    # command line.
    listing: |
      ${var qc_files_array = inputs.qc_files_array;
        var qc_files_array_of_array = inputs.qc_files_array_of_array;
        var output_array = [];

      if ( qc_files_array != null ){
        for (var i=0; i<qc_files_array.length; i++){
          output_array.push(qc_files_array[i])}}
      if ( qc_files_array_of_array != null ){
        for (var i=0; i<qc_files_array_of_array.length; i++){
          for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){
            output_array.push(qc_files_array_of_array[i][ii])}}}
      return output_array} 

baseCommand: ["multiqc"]
arguments:
  - valueFrom: --zip-data-dir
    position: 1
  - valueFrom: "'log_filesize_limit: 100000000'"
    position: 1
    prefix: --cl_config
  - valueFrom: $(runtime.outdir)
    position: 2
    prefix: --outdir
  - valueFrom: $(runtime.outdir)
    position: 4

inputs:
  qc_files_array:
    doc: |
      qc files which shall be part of the multiqc summary;
    type:
      - "null"
      - type: array
        items: File
  qc_files_array_of_array:
    doc: |
      qc files which shall be part of the multiqc summary;
    type:
      - "null"
      - type: array
        items:
          type: array
          items: File
  outdir:
    type: string?
    doc: "output directory"
    inputBinding:
      position: 1
      prefix: -o
      valueFrom: '$(runtime.outdir)'
  report_name:
    doc: name used for the html report and the corresponding zip file
    type: string?
    default: multiqc
    inputBinding:
      valueFrom: $(self + "_report")
      prefix: --filename
      position: 2
  
outputs:
  multiqc_zip:
    type: File
    outputBinding:
      glob: $(inputs.report_name + "_report_data.zip")
  multiqc_html:
    type: File
    outputBinding:
      glob: $(inputs.report_name + "_report.html")

#Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-https.rdf
  - http://edamontology.org/EDAM_1.18.owl