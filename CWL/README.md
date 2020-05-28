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
|reference|Directory where the reference genome should be put. See the ```README.md``` inside that directory|

## How to run locally

Firstly, clone this repository

```
$ git clone 
```

and then download and unzip the reference genome (see ```reference/README.md```)

All the workflows have been tested with [cwltool](https://github.com/common-workflow-language/cwltool).
[Docker](https://www.docker.com/) is also required locally.

After having installed the prerequisites and while in the Python virtual
environment created for CWL:
