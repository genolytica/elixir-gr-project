#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: Workflow

label: Complete Mapping and Quality Control pipeline for Single-end data that also executes functions of the metaseqR2 package
doc: |
  A workflow which:
  i) Runs Hisat2 on each fastq while generating a fastq file containing unmapped reads 
  ii) Runs Bowtie2 using the --very-sensitive-local option, on the unmapped reads
  iii) A subworkflow is employed to prepare BAM files, corresponding to the 
      unmapped-remaped reads, to be merged with the mapped reads form Hisat2
  iv) After the merged BAMs are created, they are sorted and assigned a user defined 
       name (samle identifier)
  v) Corresponding index file for each sample is generated.
  vi) BAMs and index files are used to create bigWig files to be used for exploring RNA signal
      in genome browsers
  vii) Above files are also used to calculate counts (counts.RData object) and generate the
       metaseqR2 html report
  viii) At the same time it performs quality control over the FASTQ using fastqc and assembles
        the MultiQC report
  
requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  MultipleInputFeatureRequirement: {}
  
inputs:
  # Input FASTQs
  fastq:
    type: File[]

  #Cores
  cores_hisat2:
    type: int?
  cores_bowtie2:
    type: int? 
  cores_samtools:
    type: int?
  cores_qc:
    type: int?
  cores_metaseqr2
    type: int?
	
  # HISAT2 alignment
  add_chr:
    type: boolean?
  hs_idx_basedir:
    type: Directory
  hs_idx_basename:
    type: string
  sam_output:
    type: string?
    default: output.sam
  unmapped_fastq:
    type: string?
    default: unmapped.fastq

  # BOWTIE2 alignment - unmapped
  local:
    type: string?
  bt_idx_basedir:
    type: Directory
  bt_idx_basename:
    type: string
  bt_output:
    type: string?
    default: bowtie.sam

  # Samtools-View arguments
  output_as_bam:
    type: boolean?
    default: true
  input_as_sam:
    type: boolean?
    default: true
  include_header:
    type: boolean?
    default: false
  bam_output:
    type: string
    default: random
  header_only:
    type: boolean?
    default: false

  # MSI - arguments
  force_overwrite:
    type: boolean?
  sorted_bam:
    type: string[] 
	
  # MetaseqR2 arguments
  targets:
    type: File  
#  targets:
#    type: 
#      type: array
#      items:
#        type: array
#        items: string
  path:
    type: Directory?  
  organism:
    type: string
  urlbase:
    type: string?
  stranded:
    type: boolean?
  normto:
    type: int?
  hubname:
    type: string?
  hubslbl:
    type: string?
  hubllbl:
    type: string?
  hubmail:
    type: string?
  exportpath:
    type: Directory?
  overwrite:
    type: boolean?
  excludelist:
    type: File?
  filetype:
    type: string?
  contrast:
    type: string
  libsizelist:
    type: File?
  idcol:
    type: int?
  gccol:
    type: int?
  namecol:
    type: int?
  btcol:
    type: int?
  annotation:
    type: string?
  refdb:
    type: string?
  version:
    type: dtring?
  translevel:
    type: string?
  counttype:
    type: string
  utrOpts_frac:
    type: float?  
  utrOpts_minlen:
    type: int?
  utrOpts_dnstrm:
    type: int?
  exonfltr:
    type: boolean?  
  exonfltr_exonsprgene
    type: int?
  exonfltr_minexons:
    type: int?
  exonfltr_frac:
    type: float?  
  genefltr:
    type: boolean?  
  genefltr1_length:
    type: int?  
  genefltr2_avgperbp:
    type: int?  
  genefltr2_avgquantile:
    type: float?  
  genefltr3_expmedian:
    type: boolean?  
  genefltr3_expmean:
    type: boolean?  
  genefltr3_expquantile:
    type: float?  
  genefltr3_expknown:
    type: File?  
  genefltr4_biotype:
    type: File?  
  genefltr5_frac:
    type: float?  
  genefltr5_mincount:
    type: int?  
  genefltr5_percon:
    type: boolean?  
  whenapplyfilter:
    type: string?
  normalization:
    type: string?
  normargs:
    type: string?
  statistics:
    type: string?
  statargs:
    type: string?
  adjmethod:
    type: string?
  metap:
    type: string?
  weight:
    type: string?
  nperm:
    type: int?
  pcut:
    type: float?
  logoffset:
    type: float?
  poffset:
    type: float?
  preset:
    type: string?
  qcplots:
    type: string?
  figformat:
    type: string?
  outlist:
    type: boolean?
  xprtwhere:
    type: Directory?
  xprtwhat:
    type: string?
  xprtscale:
    type: string?
  xprtvalues"
    type: string?
  xprtstats:
    type: string?
  xprtcountstbl:
    type: boolean?
  report:
    type: boolean?
  topreport:
    type: float?
  templatereport:
    type: string?
  genemodel:
    type: boolean?
  verbose:
    type: boolean?
  runlog:
    type: boolean?
  reportdb:
    type: string?
  localdb:
    type: Directory?
  progressfun:
    type: string?
  offlinereport:
    type: boolean?
  exportr2c:
    type: boolean?  
  
outputs:
  fastqc_zip:
    type: File[]
    outputSource: qc/fastqc_zip
  fastqc_html:
    type: File[]
    outputSource: qc/fastqc_html
  multiqc_zip:
    type: File
    outputSource: qc/multiqc_zip
  multiqc_html:
    type: File
    outputSource: qc/multiqc_html
  hisat2_sam:
    type: File[]
    outputSource: hisat2/output_sam
  hisat2_unmpd:
    type: File[]
    outputSource: hisat2/output_unmapped
  samtools_header:
    type: File[]
    outputSource: samtools_view_header/output_bam
  samtools_bam:
    type: File[]
    outputSource: samtools_view_bam/output_bam
  bowtie_sam:
    type: File[]
    outputSource: bowtie2/bowtie2_sam
  remap_bam:
    type: File[]
    outputSource: unmapped_wf/remap_bam
  final_bam:
    type: File[]
    outputSource: msi_wf/final_bam
  final_bai:
    type: File[]
    outputSource: msi_wf/final_bai
  tracks_txt:
    type: File[]
    outputSource: metaseqr2_tracks/tracks_output
  bigwig_unstranded:
    type: File[]
    outputSource: metaseqr2_tracks/bigwig
  bigwig_stranded:
    type: Directory
    outputSource: metaseqr2_tracks/bigwig2

steps:
  qc:
    run: ../workflows/qc-wf.cwl
    in:
      threads: cores_qc
      fastq: fastq
    out: [fastqc_zip,  fastqc_html, multiqc_zip, multiqc_html]

  hisat2:
    run: ../tools/hisat2_single.cwl
    scatter: fastq
    in:
      cores: cores_hisat2
      hisat2_idx_basedir: hs_idx_basedir
      hisat2_idx_basename: hs_idx_basename
      fastq: fastq
      sam_output:
        valueFrom: $(inputs.fastq.nameroot).hisat2.sam
      unmapped_fastq:
        valueFrom: $(inputs.fastq.nameroot)_unmapped.fastq
    out: [output_sam, output_unmapped]

  bowtie2:
    run: ../tools/bowtie2_single.cwl
    scatter: fq
    in:
      cores: cores_bowtie2
      local_preset: local
      bowtie2_idx_basedir: bt_idx_basedir
      bowtie2_idx_basename: bt_idx_basename
      fq: hisat2/output_unmapped
      sam_output:
        valueFrom: $(inputs.fq.nameroot)_remap_tmp.sam
    out: [bowtie2_sam]

  samtools_view_header:
    run: ../tools/samtools-view.cwl
    scatter: input_sam
    in:
      cores: cores_samtools
      header_only: header_only
      input_sam: hisat2/output_sam
      bam_output:
        valueFrom: $(inputs.input_sam.nameroot).header.sam
    out: [output_bam]

  samtools_view_bam:
    run: ../tools/samtools-view.cwl
    scatter: input_sam
    in:
      cores: cores_samtools
      output_as_bam: output_as_bam
      input_as_sam: input_as_sam
      include_header: include_header
      input_sam: hisat2/output_sam
      bam_output:
        valueFrom: $(inputs.input_sam.nameroot).bam
    out: [output_bam]

  unmapped_wf:
    run: ../workflows/unmapped.cwl
    in:
      cores: cores_samtools
      input_sam: bowtie2/bowtie2_sam
      header: samtools_view_header/output_bam
      output_as_bam: output_as_bam
      input_as_sam: input_as_sam
      include_header: include_header
    out: [remap_bam]

  msi_wf:
    run: ../workflows/merge-sort-index.cwl
    in:
      cores: cores_samtools
      force_overwrite: force_overwrite
      input_bam:
        - samtools_view_bam/output_bam
        - unmapped_wf/remap_bam
      sorted_bam: sorted_bam
    out: [final_bam, final_bai]


  metaseqr2_tracks:
    run: ../tools/reads2tracks.cwl
    in:
      targets: targets
      organism: organism
      path: path
      urlbase: urlbase
      stranded: stranded
      normto: normto
      hubname: hubname
      hubslbl: hubslbl
      hubllbl: hubllbl
      hubmail: hubmail
      exportpath: exportpath
      overwrite: overwrite
      rc: cores_metaseqr2
    out: [tracks_output, bigwig, bigwig2]
      
  metaseqr2_counts:
    run: ../tools/;redas2counts.cwl
    in:
      targets: targets
      excludelist: excludelist
      path: path
      filetype: filetype
      contrast: contrast
      libsizelist: libsizelist
      idcol: idcol
      gccol: gccol
      namecol: namecol
      btcol: btcol
      annotation: annotation
      org: organism
      refdb: refdb
      version: version
      translevel: translevel
      counttype: counttype
      utrOpts_frac  : utrOpts_frac
      utrOpts_minlen: utrOpts_minlen
      utrOpts_dnstrm: utrOpts_dnstrm
      exonfltr: exonfltr
      exonfltr_exonsprgene: exonfltr_exonsprgene
      exonfltr_minexons: exonfltr_minexons
      exonfltr_frac: exonfltr_frac
      genefltr: genefltr
      genefltr1_length: genefltr1_length
      genefltr2_avgperbp: genefltr2_avgperbp
      genefltr2_avgquantile: genefltr2_avgquantile
      genefltr3_expmedian: genefltr3_expmedian
      genefltr3_expmean: genefltr3_expmean
      genefltr3_expquantile: genefltr3_expquantile
      genefltr3_expknown: genefltr3_expknown
      genefltr4_biotype: genefltr4_biotype
      genefltr5_frac: genefltr5_frac
      genefltr5_mincount: genefltr5_mincount
      genefltr5_percon: genefltr5_percon
      whenapplyfilter: whenapplyfilter
      normalization: normalization
      normargs: normargs
      statistics: statistics
      statargs: statargs
      adjmethod: adjmethod
      metap: metap
      weight: weight
      nperm: nperm
      pcut: pcut
      logoffset: logoffset
      poffset: poffset
      preset: preset
      qcplots: qcplots
      figformat: figformat
      outlist: outlist
      xprtwhere: xprtwhere
      xprtwhat: xprtwhat
      xprtscale: xprtscale
      xprtvalues: xprtvalues
      xprtstats: xprtstats
      xprtcountstbl: xprtcountstbl
      rc: cores_metaseqr2
      report: report
      topreport: topreport
      templatereport: templatereport
      genemodel: genemodel
      verbose: verbose
      runlog: runlog
      reportdb: reportdb
      localdb: localdb
      progressfun: progressfun
      offlinereport: offlinereport
      exportr2c: exportr2c
    out: []

# Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl
