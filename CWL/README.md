# elixir-gr project-cwl-workflows

Porting of elixir-gr workflows to CWL

## Structure

The following table describes the structure of the repository.

|Directory|Contents|
|-----------|-----------|
|tools|CWL descriptions for single tools (class ```CommandLineTool```) which can be executed independently with a CWL job description (see the ```job_examples``` directory) and are used in CWL workflows (class ```Workflow```) in the ```workflows``` directory|
|workflows|CWL descriptions for SeqCVIBE workflows (class ```Workflow```) that implement various analytical steps as well as complete workflows. See the documentation (```doc``` field) inside each workflow|
|jobs|YAML file examples that provide the inputs to individual command line tools instances or complete workflows. Files are named according to the respective CWL workflow with the ```-job.yml``` suffix|
|example_data|Toy data required for testing the workflows|
|indices|Directory where HISAT2 and BOWTIE2 indices should be put. See the ```README.md``` inside that directory|

## How to run locally

Firstly, clone this repository

```
$ git clone 
```

and then download and unzip the reference genome (see ```reference/README.md```)

All the workflows have been tested with [cwltool](https://github.com/common-workflow-language/cwltool).
[Docker](https://www.docker.com/) is also required locally.

and then download and unzip the reference genome (see ```reference/README.md```)

All the workflows have been tested with [cwltool](https://github.com/common-workflow-language/cwltool).
[Docker](https://www.docker.com/) is also required locally.

After having installed the prerequisites and while in the Python virtual
environment created for CWL:

### Single tools

Make a directory to place the output 

```
$ mkdir output
```
and inside this directory, the specific directories for the metaseqR2 results

```
$ mkdir metaseqR2_results
```
```
$ mkdir signal_tracks
```

Run some tests for single tools (class ```CommandLineTool```)

1\. HISAT2 alignment for paired-end data

```
$ cwl-runner --outdir ./output \
    ./tools/hisat2_paired.cwl ./jobs/hisat2_paired-job.yml
```

2\. BOWTIE2 alignment on unmapped FASTQ output of HISAT2 

```
$ cwl-runner --outdir ./output \
    ./tools/bowtie2_paired.cwl ./jobs/bowtie2-paired-job.yml
```

3\. Convert the SAM output of HISAT2 to BAM

```
$ cwl-runner --outdir ./output \
    ./tools/samtools-view.cwl ./jobs/samtools-view-convert.yml
```

4\. Run metaseqR2 to get signal tracks (unstranded)

```
$ cwl-runner --outdir ./output \
    ./tools/reads2tracks.cwl ./jobs/r2t-job.yml
```

5\. Run metaseqR2 to get read counts and the metaseqR2 html report

```
$ cwl-runner --outdir ./output \
     ./tools/reads2counts.cwl ./jobs/r2c-job.yml   
```

### Workflows

Run some tests for workflows (class ```Workflow```)

1\. Samtools post-alignment processing on BOWTIE2 output for preparing BAM to be merged with HISAT2 alignments

```
$ cwl-runner --outdir ./output \
    ./workflows/unmapped.cwl ./jobs/unmapped-job.yml
```

2\. Mapping & Quality Control workflow

```
$ cwl-runner --outdir ./output \
    ./workflows/mapping-pe-qc.cwl ./jobs/mapping-pe-qc-job.yml
```

3\. Complete Mapping, Quality Control and Quantification workflow

```
$ cwl-runner --outdir ./output \
    ./workflows/mapping-pe-qc-r.cwl ./jobs/mapping-pe-qc-r-job.yml
```