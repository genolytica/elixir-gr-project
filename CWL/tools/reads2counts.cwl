#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: metaseqR2 for acquiring Counts
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
      path: ../runMetaseqR2Counts.R
    inputBinding:
      position: 1
  targets:
#    type: Any
#    type:
#      type: array
#      items: string[]
#    type: 
#      type: array
#      items:
#        type: array
#        items: string
    type: File
    doc: "targets file equivalent"
    inputBinding:
      prefix: --samplelist
      position: 2
#  targets:
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
  excludelist:
    type: File?
    doc: "A list of samples to exclude, in the same format as the input targets.txt"
    inputBinding:
     prefix: --excludelist
     position: 3
  filetype:
    type: string?
    doc: "The type of raw input files"
    inputBinding:
      prefix: --filetype
      position: 4
  path:
    type: Directory?
    doc: "Path where all the BED/BAM files are placed"
    default:
      class: Directory
      path: ../example_data
    inputBinding:
      prefix: --path
      position: 5
  contrast:
    type: string
    doc: "Example: condition1_vs_condition2. Condition1 MUST be the control"
    inputBinding:
      prefix: --contrast
      position: 6
  libsizelist:
   type: File?
   doc: "An optional named list where names represent samples and members are the library sizes for each sample"
   inputBinding:
     prefix: --libsizelist
     position: 7
  idcol:
    type: int?
    doc: "Column number where the unique gene accessions are. Default to 4"
    inputBinding:
      prefix: --idcol
      position: 8
  gccol: 
    type: int?
    doc: "Column number where each gene's GC content is given. If not provided, GC content normalization provided by EDASeq will not be available"
    inputBinding:
      prefix: --gccol
      position: 9
  namecol:
      type: int?
      doc: "Column number where the HUGO gene symbols are given. If not provided, it will not be available when reporting results. In addition, the 'known' gene filter will not be available for application"
      inputBinding:
        prefix: --namecol
        position: 10
  btcol:
    type: int?
    doc: "Column number where the gene biotypes are given. If not provided, the 'biodetection', 'countsbio', 'saturation', 'filtered' and 'biodist' plots will not be available"
    inputBinding:
      prefix: --btcol
      position: 11
  annotation:
    type: string?
    doc: "NULL (default), embedded or a list with a path to a GTF file and certain required metadata"
    inputBinding:
      prefix: --annotation
      position: 12
  organism:
    type: string
    doc: "organism annotation"
    inputBinding:
      prefix: --org
      position: 13
  refdb:
    type: string?
    default: ensembl
    doc: "The reference annotation repository from which to retrieve annotation elements: ensembl(default), ucsc or refseq"
    inputBinding:
      prefix: --refdb
      position: 14
  version:
    type: string?
    doc: "An integer denoting the version of the annotation to use"
    inputBinding:
      prefix: --version
      position: 15
  translevel:
    type: string?
    doc: "Transcriptional unit for DE : gene (default) or transcript"
    inputBinding:
      prefix: --translevel
      position: 16
  counttype:
    type: string
    doc: "The type of reads inside the counts file : gene, exon or utr"
    inputBinding:
      prefix: --counttype
      position: 17
  utropts_frac:
    type: float?
    default: 1
    doc: "The fraction (0-1) of the 3' UTR region to count reads in. Default 1"
    inputBinding:
      prefix: --utrOpts_frac
      position: 18
  utropts_minlen:
    type: int?
    default: 300
    doc: "The minimum acceptable 3'UTR length irrespective of utrOpts_frac argument. Default 300"
    inputBinding:
      prefix: --utrOpts_minlen
      position: 19
  utropts_dnstrm:
    type: int?
    default: 50
    doc: "The number of base pairs to flank the end of the 3' UTR of transcripts. Default 50"
    inputBinding:
      prefix: --utrOpts_dnstrm
      position: 20
  exonfltr:
    type: boolean?
#    default: true
    doc: "Set exonfltr FALSE to NOT apply the minActiveExons filter"
    inputBinding:
      prefix: --exonfltr
      position: 21
  exonfltr_exonsprgene:
    type: int?
    default: 5
    doc: "minActiveExons filter: Exons per gene. Defaults to 5"
    inputBinding:
      prefix: --exonfltr_exonsprgene
      position: 22
  exonfltr_minexons:
    type: int?
    default: 2
    doc: "minActiveExons filter: read presence is required in at least exonfltr_minexons of exonfltr_exonsprgene. Defaults to 2"
    inputBinding:
      prefix: --exonfltr_minexons
      position: 23
  exonfltr_frac:
    type: float?
    default: 0.2
    doc: "Read presence is required in a exonfltr_frac fraction of the total exons. Defaults to 1/5"
    inputBinding:
      prefix: --exonfltr_frac
      position: 24
  genefltr:
    type: boolean?
#    default: true
    doc: "Set genefltr=FALSE to NOT apply any gene filtering"
    inputBinding:
      prefix: --genefltr
      position: 25
  genefltr1_length:
    type: int?
    default: 500
    doc: "Gene filter1: length filter where genes are accepted for further analysis if they are above genefltr_length. Defaults to 500"
    inputBinding:
      prefix: --genefltr1_length
      position: 26
  genefltr2_avgperbp:
    type: int?
    default: 100
    doc: "Gene filter2: a gene is accepted for further analysis if it has more average reads than the genefltr_avgquantile of the average count distribution per genefltr_avgperbp base pairs. Defaults to 100"
    inputBinding:
      prefix: --genefltr2_avgperbp
      position: 27
  genefltr2_avgquantile:
    type: float?
    default: 0.25
    doc: "Gene filter2: a gene is accepted for further analysis if it has more average reads than the genefltr_avgquantile of the average count distribution per genefltr_avgperbp base pairs. Defaults to 0.25"
    inputBinding:
      prefix: --genefltr2_avgquantile
      position: 28
  genefltr3_expmedian:
    type: boolean?
#    default: true
    doc: "Gene filter3: based on the overall expression of a gene. Genes below the median of the overall count distribution are not accepted for further analysis. Defaults to TRUE"
    inputBinding:
      prefix: --genefltr3_expmedian
      position: 29
  genefltr3_expmean:
    type: boolean?
#    default: false
    doc: "Gene filter3: based on the overall expression of a gene. Genes below the mean of the overall count distribution are not accepted for further analysis. Defaults to FALSE"
    inputBinding:
      prefix: --genefltr3_expmean
      position: 30
  genefltr3_expquantile:
    type: float?
    doc: "Gene filter3: based on the overall expression of a gene. Genes below the specified quantile of the total counts distribution are not accepted for further analysis"
    inputBinding:
      prefix: --genefltr3_expquantile
      position: 31
  genefltr3_expknown:
    type: File?
    doc: "Gene filter3: based on the overall expression of a gene. A set of known not-expressed genes in the system under investigation are used to estimate an expression cutoff"
    inputBinding:
      prefix: --genefltr3_expknown
      position: 32
  genefltr4_biotype:
    type: File?
    doc: "Gene filter4: genes with a certain biotype (MUST be contained in the annotation,thus it's better to use annotation='download') are excluded from the analysis"
    inputBinding:
      prefix: --genefltr4_biotype
      position: 33
  genefltr5_frac:
    type: float?
    default: 0.25
    doc: "Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to 0.25"
    inputBinding:
      prefix: --genefltr5_frac
      position: 34
  genefltr5_mincount:
    type: int?
    default: 10
    doc: "Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to 10"
    inputBinding:
      prefix: --genefltr5_mincount
      position: 35
  genefltr5_percon:
    type: boolean?
#    default: false
    doc: "Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to FALSE"
    inputBinding:
      prefix: --genefltr5_percon
      position: 36
  whenapplyfilter:
    type: string?
    default: postnorm
    doc: "When to apply exon and/or gene filtering: postnorm (default) or prenorm"
    inputBinding:
      prefix: --whenapplyfilter
      position: 37
  normalization:
    type: string?
    default: deseq
    doc: "The normalization algorithm to be applied on the count data"
    inputBinding: 
      prefix: --normalization
      position: 38
  normargs:
    type: string?
#    default: NULL
    doc: "A named list whose names are the names of the normalization algorithm parameters and its members parameter values. Leave NULL for the defaults of normalization"
    inputBinding:
      prefix: --normargs
      position: 39
  statistics:
    type: string?
    default: deseq
    doc: "One or more statistical analyses to be performed by the metaseqr2 pipeline. Enter as comma-separated list. Default is deseq"
    inputBinding:
      prefix: --statistics
      position: 40 
  statargs:
    type: string?
    doc: "A named list whose names are the names of the statistical algorithms used in the pipeline. Leave NULL for the defaults of statistics"
    inputBinding:
      prefix: --statargs
      position: 41
  adjmethod:
    type: string?
    default: BH
    doc: "The multiple testing p-value adjustment method. It can be one of p.adjust.methods or 'qvalue' from the qvalue Bioconductor package. Defaults to 'BH' for Benjamini-Hochberg correction"
    inputBinding:
      prefix: --adjmethod
      position: 42
  metap:
    type: string?
    default: simes
    doc: "The meta-analysis method to combine p-values from multiple statistical tests. It can be one of 'simes'(default), 'bonferroni', 'minp', 'maxp', 'weight, 'pandora', 'dperm_min', 'dperm_max', 'dperm_weight', 'fisher', 'fperm', 'whitlock' or 'none'"
    inputBinding:
      prefix: --metap
      position: 43
  weight:
    type: string?
    doc: "A vector of weights with the same length as the statistics vector containing a weight for each statistical test. It should sum to 1"
    inputBinding:
      prefix: --weight
      position: 44
  nperm: 
    type: int?
    default: 10000
    doc: ""
    inputBinding:
      prefix: --nperm
      position: 45
  pcut:
    type: float?
    doc: "A p-value cutoff for exporting differentially genes, default is to export all the non-filtered genes"
    inputBinding:
      prefix: --pcut
      position: 46
  logoffset:
    type: float?
    default: 1
    doc: "An offset to be added to values during logarithmic transformations in order to avoid Infinity (default is 1)"
    inputBinding:
      prefix: --logoffset
      position: 47
  poffset:
    type: float?
    doc: "A value between 0 and 1 to multiply potential zero p-values with for the combination methods including weighting or NULL (default)"
    inputBinding:
      prefix: --poffset
      position: 48
  preset:
    type: string?
    doc: "One of 'all_basic', all_normal', 'all_full', 'medium_basic', 'medium_normal', 'medium_full', 'strict_basic', 'strict_normal' or 'strict_full', each of which control. The strictness of the analysis and the amount of data to be exported"
    inputBinding:
      prefix: --preset
      position: 49
  qcplots:
    type: string?
    default: mds,biodetection,rnacomp,countsbio,saturation,readnoise,filtered,boxplot,gcbias,lengthbias,meandiff,meanvar,deheatmap,volcano,mastat,biodist,statvenn,foldvenn
    doc: "A set of diagnostic plots to show/create"
    inputBinding:
      prefix: --qcplots
      position: 50
  figformat:
    type: string
    default: png,jpg,tiff,bmp,pdf,ps
    doc: "The format of the output diagnostic plots: png,jpg,tiff,bmp,pdf,ps"
    inputBinding:
      prefix: --figformat
      position: 51
  outlist:
    type: boolean?
#   default: false
    doc: "Export a list with the results? Default FALSE"
    inputBinding:
      prefix: --outlist
      position: 52 
  xprtwhere:
    type: Directory?
    doc: "An output directory for the project results"
#    default:
#      class: Directory
#      path: ../output
    inputBinding:
      prefix: --xprtwhere
      position: 53
  xprtwhat:
    type: string?
    default: annotation,p_value,adj_p_value,meta_p_value,adj_meta_p_value,fold_change,stats,counts,flags
    doc: "The content of the final lists"
    inputBinding:
      prefix: --xprtwhat
      position: 54
  xprtscale:
    type: string?
    default: natural,log2,log10,vst,rpgm
    doc: "Export values from one or more transformations applied to the data. Enter as comma-separated list to select 1 or more"
    inputBinding:
      prefix: --xprtscale
      position: 55
  xprtvalues:
    type: string?
    default: raw,normalized
    doc: "It can be one or more of 'raw' to export raw values (counts etc.) and 'normalized' to export normalized counts. Enter as comma-separated list"
    inputBinding:
      prefix: --xprtvalues
      position: 56
  xprtstats:
    type: string?
    default: mean,median,sd,mad,cv,rcv
    doc: "Calculate and export several statistics on raw and normalized counts, condition-wise. Enter as comma-separated list to select 1 or more"
    inputBinding:
      prefix: --xprtstats
      position: 57
  xprtcountstbl:
    type: boolean?
#    default: false
    doc: "Exports the calculated read counts table when input is read from bam files and the normalized count table in all cases. Defaults to FALSE"
    inputBinding:
      prefix: --xprtcountstbl
      position: 58
  rc:
    type: float?
    doc: "Fraction of cores to use"
    inputBinding:
      prefix: --rc
      position: 59
  report:
    type: boolean?
#    default: true
    doc: "A logical value controlling whether to produce a summary report or not. Defaults to TRUE"
    inputBinding:
      prefix: --report
      position: 60
  toreport:
    type: float?
    default: 0.1
    doc: "A fraction of top statistically significant genes to append to the HTML report"
    inputBinding:
      prefix: --topreport
      position: 61
  templatereport:
    type: string?
    doc: "An HTML template to use for the report. Do not change this unless you know what you are doing"
    inputBinding:
      prefix: --templatereport
      position: 62
  genemodel:
    type: boolean?
#    default: true
    doc: "In case of exon analysis, a list with exon counts for each gene will be saved to the file xprtwhere/data/gene_model.RData"
    inputBinding:
      prefix: --genemodel
      position: 63
  verbose:
    type: boolean?
#    default: true
    doc: "Print informative messages during execution? Defaults to TRUE"
    inputBinding:
      prefix: --verbose
      position: 64
  runlog:
    type: boolean?
#    default: true
    doc: "Write a log file of the metaseqr2 run using package log4r. Defaults to TRUE"
    inputBinding:
      prefix: --runlog
      position: 65
  reportdb:
    type: string?
    default: sqlite
    doc: "Database system to use for storing the report intereactive graphs. Can be 'sqlite' (default) or 'dexie'"
    inputBinding:
      prefix: --reportdb
      position: 66
  localdb:
    type: Directory?
    doc: "The metaseqR2 annotation database location"
    inputBinding:
      prefix: --localdb
      position: 67
  progressfun:
    type: string?
    doc: "A function which updates a Progress object from shiny. This function must accept a detail argument. See http://shiny.rstudio.com/articles/progress.html"
    inputBinding:
      prefix: --progressfun
      position: 68
  offlinereport:
    type: boolean?
#    default: true
    doc: "TRUE (default) to download and include the required JavaScript libraries to properly view the report offline. Ignored if report=FALSE"
    inputBinding:
      prefix: --offlinereport
      position: 69
  exportr2c:
    type: boolean?
#    default: true
    doc: "Exports the counts.rda object. Defaults to TRUE"
    inputBinding:
      prefix: --exportr2c
      position: 70

outputs:
  data:
    type: File[]
    outputBinding:
      glob: "*.RData"
  javascript:
    type: File[]
    outputBinding:
      glob: "*.js"
  results:
    type: Directory
    outputBinding:
      glob: "*result*"
  results2:
    type: Directory?
    outputBinding:
      glob: $(inputs.xprtwhere)

#Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf
  - http://edamontology.org/EDAM_1.18.owl