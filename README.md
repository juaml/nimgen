# NImGen

NImGen is a Python-based extensible open source toolbox for elucidating
neuroimaging-genetics mappings. The toolbox allows extracting gene expression
levels for parcellations from the [Allen Human Brain Atlas (AHBA)](https://portal.brain-map.org/) with the help
of the [abagen](https://abagen.readthedocs.io/en/stable/) package. 
Each gene’s expression levels can be correlated to an MRI
marker of interest. With the assistance of a gene enrichment analysis
tool [(WebGestalt)](http://www.webgestalt.org/), the significant genes are examined to obtain higher-order
biological processes.  To properly account for spatial autocorrelation in the brain,
and because increasing the number of parcels can inflate statistical significance,
NImGen uses a permutation-testing approach based on the [BrainSMASH](https://brainsmash.readthedocs.io/en/latest/)
and [neuromaps](https://netneurolab.github.io/neuromaps/) python packages.
An efficient, reproducible, and rapid pipeline with several MRI markers may be
constructed due to the package's modular design and functionality.


## Installing

NImGen was developed and tested using Python 3.9. Backwards-compatibility with other Python3 versions is expected but not guaranteed.  
NImGen is most easily installed via GitHub repository  

In order for the nimgen pipeline to interface with [the WebGestalt R package](https://cran.r-project.org/web/packages/WebGestaltR/index.html),
you will of course need R and the WebGestalt dependency. Now, if you do not have these,
don't worry. The pipeline will still run as is and give the desired output, except you will not
get output for gene enrichment analysis. However, you can use the output to
perform a gene enrichment analysis yourself, if you wish to do so of cource.

If you do want nimgen to perform the gene enrichment analysis via WebGestalt
for you, then you may want to use a [conda environment](https://docs.conda.io/en/latest/miniconda.html)
for your installation.

For example, to create and activate a conda environment with the WebGestalt R dependencies, run:

```
conda create -n my_nimgen_environment r-essentials r-base r-WebGestaltR python=3.9
conda activate my_nimgen_environment
```

Either way, to install the nimgen package, you can simply install from GitHub using pip:

```
pip install git+https://github.com/juaml/nimgen.git
```

## Usage

Currently, the main pipeline provided by nimgen is intended to run on an [HTCondor-based cluster](https://htcondor.org/),
but support for parallelisation on local computers or other platforms may still follow.
To run the nimgen pipeline, you need two input files and one configuration file.
In this quick tutorial, we will show you how to use a basic nimgen pipeline.

### Inputs

1. Marker: A brain image/map as NIfTI (.nii) file in MNI152 space. An example could be voxel based morphometry as e.g. outputted by CAT.
2. Atlas/Parcellation: Brain atlases i.e. Schaefer atlas 100, 500 or 1000 as NIfTI (.nii), again in MNI152 space.  

### Configuration

The pipeline is configured via a `.yaml` file.
A basic configuration file for nimgen may look as follows:

```
name: example_pipeline
verbosity: INFO

seed: 100
pipeline:
    type: "HTCondor"
    step_1:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
    step_2:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
    step_3:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
    step_4:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB

conda_env: "my_nimgen_environment"

markers:
    - path: "example_marker.nii"
      parcellation: 
        - "example_atlas.nii"

n_surrogate_maps: 100

correlation_method:
    - spearman

alpha:
    - 0.05
```

Let's quickly discuss each field and its signifiance. At the top part of the file you see 
these three fields:

```
name: example_pipeline
verbosity: INFO

seed: 100
```

The `name` field specifies the name of your pipeline and can be any arbitrary string that you like.
In this case, the configuration file will tell nimgen to create a folder called `example_pipeline` with
all pipeline related files inside it.
The `verbosity` field will tell nimgen how verbose logging output should be when running the individual pipeline
steps. Possible values here are any of the following: `{DEBUG,INFO,WARNING,ERROR}`.
The `seed` field sets a random seed and is important for the creation of permutation-based null maps,
as of course, randomness is involved in this step.

Next, you can see a rather large block, specifying some more details of the pipeline.
Importantly, you can see the `type` key. Currently, this can only be set to "HTCondor",
but other pipeline types will hopefully soon follow.
Below, you can see `CPU`, `MEMORY`, and `DISK` specifications for each pipeline step.
Their values should be a number with a unit as displayed below.
We will talk about the pipeline steps shortly.

```
pipeline:
    type: "HTCondor"
    step_1:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
    step_2:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
    step_3:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
    step_4:
        CPU: 1
        MEMORY: 8GB
        DISK: 10GB
```

Below this larger pipeline block, you can see a field specifying a conda environment:

```
conda_env: "my_nimgen_environment"
```

However, if you dont use conda, you can also specify a path to a python venv as follows:
(Note that `path_to_venv` is the field/key and `path/to/my/project/venv` is the actual path to source the venv)

```
path_to_venv: path/to/my/project/venv
```

In addition, you may also explicitly want to set your path to `Rscript` for the `WebGestaltR` gene enrichment analysis.
By default this will be set to `/usr/bin/env Rscript`.

Lastly, you can see this block in the yaml file:

```
markers:
    - path: "example_marker.nii"
      parcellation: 
        - "example_atlas.nii"

n_surrogate_maps: 100

correlation_method:
    - spearman

alpha:
    - 0.05
```

The `markers` field specifies a list of markers (note the dash before the `path` key indicating a list.)
Each marker in turn is a dictionary, consisting of a path to the marker (this is the `path` key) and a 
list of parcellations. That is, for each marker, the pipeline will generate output for every parcellation listed there.
The list of parcellations is simply a list of paths to each parcellation. Note that each path is evaluated in the working
directory in which you run `nimgen create` (we will talk about this command soon). 

The field `n_surrogate_maps` specifies how many null maps you want to create for empirical p-value computation.
The `correlation_method` field takes a list of correlation methods and will output results for each.
Currently supported correlation methods are: `{spearman, pearson}`.

Lastly, the `alpha` field gives a list of alpha levels of significance for which gene enrichment analysis can be run.

### Run a basic pipeline

As mentioned above, to run a nimgen pipeline you need at least two inputs: An atlas file and a marker file (both are NIfTI files in MNI152 space.)
Now, say you have two such files, and save them with your above yaml configuration in your current working directory, then
running `ls` will yield the following output:

```
example_atlas.nii  example_config.yaml  example_marker.nii
```

In order to prepare a pipeline simply run `nimgen create example_config.yaml`.
When running `ls`, you will see that nimgen created an additional folder in your
current working directory:

`example_atlas.nii  example_config.yaml  example_marker.nii  example_pipeline `

If you run `tree example_pipeline`, you can see the contents of that directory:

```bash
example_pipeline
├── all_gene_outputs
├── example_config.yaml
├── markers
│   └── example_marker
│       ├── example_marker.nii
│       ├── nullmaps
│       │   └── example_atlas
│       │       └── nullmaps_results
│       └── outputs
│           └── example_atlas
├── parcellations
│   └── example_atlas
│       └── example_atlas.nii
└── submit_files
    ├── logs
    ├── nimgen.dag
    ├── run_in_venv.sh
    ├── step_1.submit
    ├── step_2.submit
    ├── step_3.submit
    └── step_4.submit

12 directories, 9 files
```

The `all_gene_outputs` folder is used to collect all output specific to all
genes used in the AHBA, and as such are not specific to any one marker.
The `markers` directory collects intermediate data such as `nullmaps` as well
as marker specific output for each marker specified in the configuration file.
Importantly, the `submit_files` directory collects all condor_submit files and
and executable called `run_in_venv.sh`. This executable is used to ensure that
each job is executed in the specified environment. To run all jobs, you can
simply go ahead and submit the `nimgen.dag` file:
```
cd example_pipeline/submit_files
condor_submit_dag nimgen.dag
```

Alternatively, you can also submit each individual step yourself. For example:
```
condor_submit step1.submit
```
Or to submit only one set of parameters for each step, have a peak into each submit file and
look at the commands. For example from `step1.submit` you can extract the following command
and run it (make sure you are in the correct environment):

```
nimgen -v INFO step_1 ../../example_pipeline/parcellations/example_atlas/example_atlas.nii
```

## Core development team
Kaustubh Patil, Fede Raimondo, Leonard Sasse, Talip Yasir Demirtaş

## Support
If you run into a problem or identify a bug in the source code, or if you have any
questions about nimgen's internals or its output, 
please send the authors an email or create a new issue on GitHub.

## References

* [Generative modeling of brain maps with spatial autocorrelation](https://www.sciencedirect.com/science/article/pii/S1053811920305243)
* [Standardizing workflows in imaging transcriptomics with the abagen toolbox](https://www.biorxiv.org/content/10.1101/2021.07.08.451635v1)
* [A practical guide to linking brain-wide gene expression and neuroimaging data](https://www.sciencedirect.com/science/article/abs/pii/S1053811919300114?via%3Dihub)
* [An anatomically comprehensive atlas of the adult human brain transcriptome](https://www.nature.com/articles/nature11405)
* [neuromaps: structural and functional interpretation of brain maps](https://www.nature.com/articles/s41592-022-01625-w)
* [WebGestalt 2017: a more comprehensive, powerful, flexible and interactive gene set enrichment analysis toolkit](https://academic.oup.com/nar/article/45/W1/W130/3791209?login=false)
