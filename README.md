# NImGen

## Requirements:

* pandas
* abagen
* selenium (only for webgestalt)

## Installing

```
git clone https://github.com/juaml/nimgen.git
cd nimgen
python setup.py develop
```

## Allen Brain Atlas 

Please get the `normalized microarray data` for all the subjects from here:  https://human.brain-map.org/static/download

Put all of them in a single directory and decompress. The directory should look like as following consisting of the six subjects:

```
├── normalized_microarray_donor10021
│   ├── MicroarrayExpression.csv
│   ├── Ontology.csv
│   ├── PACall.csv
│   ├── Probes.csv
│   ├── Readme.txt
│   └── SampleAnnot.csv
├── normalized_microarray_donor12876
│   ├── MicroarrayExpression.csv
│   ├── Ontology.csv
│   ├── PACall.csv
│   ├── Probes.csv
│   ├── Readme.txt
│   └── SampleAnnot.csv
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
