#!/usr/bin/ cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: metaseqR2 for acquiring Counts
doc: |

requirements:
  DockerRequirement:
    dockerPull: "daphnelettos/metaseqr2:1.1"

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
    type: File
    doc: "targets file"
    inputBinding:
      prefix: --samplelist
      position: 2
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
    type: Directory
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
  annotation:
    type: string?
    doc: "NULL (default), embedded or a list with a path to a GTF file and certain required metadata"
    inputBinding:
      prefix: --annotation
      position: 8
  organism:
    type: string
    doc: "organism annotation"
    inputBinding:
      prefix: --org
      position: 9
#  refdb:
#  type: string?
#  default: ensembl
#  doc: "The reference annotation repository from which to retrieve annotation elements: ensembl(default), ucsc or refseq"
#  inputBinding:
#    prefix: --refdb
#    position: 10
#	make_option(
# 		opt_str="--version",
#		action="store",
#		type="numeric",
#		default="auto",
#		help=paste0(
#			"An integer denoting the version of the annotation to use from the local annotation\n",
#			"database or fetch on the fly. For Ensembl, it corresponds to Ensembl releases, while\n",
#			"for UCSC/RefSeq, it is the date of creation (locally)."
  translevel:
    type: string?
    doc: "Transcriptional unit for DE : gene (default) or transcript"
    inputBinding:
      prefix: --translevel
      position: 11
  counttype:
    type: string
    doc: "The type of reads inside the counts file : gene, exon or utr"
    inputBinding:
      prefix: --counttype
      position: 12
  utropts_frac:
    type: int?
    default: 1
    doc: "The fraction (0-1) of the 3' UTR region to count reads in. Default 1"
    inputBinding:
      prefix: --utrOpts_frac
      position: 13
  utropts_minlen:
    type: int?
    default: 300
    doc: "The minimum acceptable 3'UTR length irrespective of utrOpts_frac argument. Default 300"
    inputBinding:
      prefix: --utrOpts_minlen
      position: 14
  utropts_dnstrm:
    type: int?
    default: 50
    doc: "The number of base pairs to flank the end of the 3' UTR of transcripts. Default 50"
    inputBinding:
      prefix: --utrOpts_dnstrm
      position: 15
  exonfltr:
    type: boolean?
    default: true
    doc: "Set exonfltr FALSE to NOT apply the minActiveExons filter"
    inputBinding:
      prefix: --exonfltr
      position: 16
  exonfltr_exonsprgene:
    type: int?
    default: 5
    doc: "minActiveExons filter: Exons per gene. Defaults to 5"
    inputBinding:
      prefix: --exonfltr_exonsprgene
      position: 16
  exonfltr_minexons:
    type: int?
    default: 2
    doc: "minActiveExons filter: read presence is required in at least exonfltr_minexons of exonfltr_exonsprgene. Defaults to 2"
    inputBinding:
      prefix: --exonfltr_minexons
      position: 17
#  exonfltr_frac:
#    type: double?
#    default: 1/5
#    doc: "Read presence is required in a exonfltr_frac fraction of the total exons. Defaults to 1/5"
#    inputBinding:
#      prefix: --exonfltr_frac
#      position: 18
  genefltr:
    type: boolean?
    default: true
    doc: "Set genefltr=FALSE to NOT apply any gene filtering"
    inputBinding:
      prefix: --genefltr
      position: 19
  genefltr1_length:
    type: int?
    default: 500
    doc: "Gene filter1: length filter where genes are accepted for further analysis if they are above genefltr_length. Defaults to 500"
    inputBinding:
      prefix: --genefltr1_length
      position: 20
  genefltr2_avgperbp:
    type: int?
    default: 100
    doc: "Gene filter2: a gene is accepted for further analysis if it has more average reads than the genefltr_avgquantile of the average count distribution per genefltr_avgperbp base pairs. Defaults to 100"
    inputBinding:
      prefix: --genefltr2_avgperbp
      position: 21
  genefltr2_avgquantile:
    type: double?
    default: 0.25
    doc: "Gene filter2: a gene is accepted for further analysis if it has more average reads than the genefltr_avgquantile of the average count distribution per genefltr_avgperbp base pairs. Defaults to 0.25"
    inputBinding:
      prefix: --genefltr2_avgquantile
      position: 22
  genefltr3_expmedian:
    type: boolean?
    default: true
    doc: "Gene filter3: based on the overall expression of a gene. Genes below the median of the overall count distribution are not accepted for further analysis. Defaults to TRUE"
    inputBinding:
      prefix: --genefltr3_expmedian
      position: 23
  genefltr3_expmean:
    type: boolean?
    default: false
    doc: "Gene filter3: based on the overall expression of a gene. Genes below the mean of the overall count distribution are not accepted for further analysis. Defaults to FALSE"
    inputBinding:
      prefix: --genefltr3_expmean
      position: 24
#  genefltr3_expquantile:
#    type: 
#		opt_str="--genefltr3_expquantile",
#		action="store",
#		default=NA,
#		help=paste0(
#			"Gene filter3: based on the overall expression of a gene. Genes below the specified\n",
#			"quantile of the total counts distribution are not accepted for further analysis."
#  genefltr3_expknown:
#    type: 
#	make_option(
#		opt_str="--genefltr3_expknown",
#		action="store",
#		default=NA,
#		help=paste0(
#			"Gene filter3: based on the overall expression of a gene. A set of known not-expressed\n",
#			"genes in the system under investigation are used to estimate an expression cutoff.\n",
#			"The value of this filter is a character vector of HUGO gene symbols (MUST be contained\n",
#			"in the annotation, thus it's better to use annotation='download') whose counts are used\n",
#			"to build a 'null' expression distribution. The 90th quantile of this distribution is then."
#		opt_str="--genefltr4_biotype",
#		action="store",
#		default="getDefaults('biotype.filter',org[1])",
#		help=paste0(
#			"Gene filter4: genes with a certain biotype (MUST be contained in the annotation,\n",
#			"thus it's better to use annotation='download') are excluded from the analysis.\n",
#			"This filter is a named list of logical, where names are the biotypes in each genome\n",
#			"and values are TRUE or FALSE. If the biotype should be excluded, the value should be TRUE else FALSE."
#		)
#	),
  genefltr5_frac:
    type: double?
    default: 0.25
    doc: "Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to 0.25"
    inputBinding:
      prefix: --genefltr5_frac
      position: 25
  genefltr5_mincount:
    type: int?
    default: 10
    doc: "Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to 10"
    inputBinding:
      prefix: --genefltr5_mincount
      position: 26 
  genefltr5_percon:
    type: boolean?
    default: false
    doc: "Gene filter5: a gene is further considered for statistical testing if genefltr5_frac (x100 for a percentage value) have more than genefltr5_mincount reads across all samples (genefltr5_percon=FALSE) or across the samples of each condition (genefltr5_percon=TRUE). Defaults to FALSE"
    inputBinding:
      prefix: --genefltr5_percon
      position: 27
  whenapplyfilter:
    type: string?
    default: postnorm
    doc: "When to apply exon and/or gene filtering: postnorm (default) or prenorm"
    inputBinding:
      prefix: --whenapplyfilter
      position: 28
  normalization:
    type: string?
    default: deaeq
    doc: "The normalization algorithm to be applied on the count data"
    inputBinding: 
      prefix: --normalization
      position: 29
  statistics:
    type: string?
    default: deseq
    doc: "Statistical analysis to be performed"
    inputBinding:
      prefix: --statistics
      position: 30 
  qcplots:
    type: string?
    doc: "mds, biodetection, rnacomp, countsbio, saturation, readnoise, filtered, boxplot, gcbias, lengthbias, meandiff, meanvar, deheatmap, volcano, mastat, biodist, statvenn, foldvenn. Defaults to qcPlots=NULL"
    inputBinding:
      prefix: --qcplots
      position: 31
  figformat:
    type: string
    doc: "The format of the output diagnostic plots: png,jpg,tiff,bmp,pdf,ps"
    inputBinding:
      prefix: --figformat
      position: 32
  outlist:
    type: boolean?
    default: false
    doc: "Export a list with the results? Default FALSE"
    inputBinding:
      prefix: --outlist
      position: 33 
  xprtwhere:
    type: Directory?
    doc: "An output directory for the project results"
    default:
      class: Directory
      path: ../output
    inputBinding:
      prefix: --xprtwhere
      position: 34
  rc:
    type: int?
    doc: "Fraction of cores to use"
    inputBinding:
      prefix: --rc
      position: 35


outputs:
  data:
    type: Directory?
    outputBinding:
      glob: "data"
#  plots:
#    type: Directory?
#    outputBinding:
#      glob: "plots"



#Metadata
$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf
  - http://edamontology.org/EDAM_1.18.owl