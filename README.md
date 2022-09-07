# NImGen

## Requirements:

* abagen
* brainsmash
* nibabel
* nilearn
* pandas
* numpy
* pingouin
* scikit_learn
* scipy
* statsmodels


## Installing

```
git clone https://github.com/juaml/nimgen.git
cd nimgen
python setup.py develop
```

## Usage

Three components are needed to use the functionality:

```
# a 3D nifti file containing ROIs labelled as integers
atlas = './data/Power_5mm.nii'
# a plain text file containing weights for each of the ROIs (ordered according to the ROI integer in the atlas)
weights = pd.read_csv('./data/rehoavg.txt', header=None)
# a directory containing the microarray data (see below)
allen_data_dir = './data/allen'
```

Once these are ready you can get the multiple-testing corrected (default BH FDR) P-values and genes as follows:

```
all_genes, sign_genes = nimgen.get_gene_expression(
    weights, atlas, allen_data_dir=allen_data_dir, save_expressions=True)
print(sign_genes)
```

The `save_expressions=True` saves the extracted expression values so that they can be reused. Please check the respective functions for more details.

The `examples` directory contains more usage examples.

Check the `abagen` library for more details: https://github.com/rmarkello/abagen
  

## Allen Brain Atlas 

Please get the `normalized microarray data` for all the subjects from here:  https://human.brain-map.org/static/download

Put all of them in a single directory and decompress. The directory should have the following structure, consisting of the six subjects and six files per subjects:

```
├── normalized_microarray_donor10021
│   ├── MicroarrayExpression.csv
│   ├── Ontology.csv
│   ├── PACall.csv
│   ├── Probes.csv
│   ├── Readme.txt
│   └── SampleAnnot.csv
├── normalized_microarray_donor12876
│   ├── ...
├── normalized_microarray_donor14380
│   ├── ...
├── normalized_microarray_donor15496
│   ├── ...
├── normalized_microarray_donor15697
│   ├── ...
└── normalized_microarray_donor9861
    ├── ...
```

You will have to provide this directory to the functions calls that get the gene expression data.

## Using webgestalt

To use webgestalt, the respective selenium driver must be installed.

Firefox drivers: https://sites.google.com/a/chromium.org/chromedriver/downloads
Firefox: https://github.com/mozilla/geckodriver/releases

Just download the respective file, unzip and place in the users PATH (`/bin` or `/usr/bin`)

Note for Macos users: 
Open a Terminal, copy the driver to `/usr/local/bin` and then execute. An error message will appear. Press close or cancel. Open system preferences, security and press the button to allow for the execution. Go to the terminal, execute again. Another warning will appear. Allow.

