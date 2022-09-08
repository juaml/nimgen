# NImGen

NImGen is a Python-based extensible open source toolbox for elucidating neuroimaging-genetics mappings. The toolbox allows extracting gene expression levels for parcellations from the Allen Human Brain Atlas (AHBA) with the help of the abagen package. Each gene’s expression levels can be correlated to an MRI marker of interest. Besides correlation analysis, the NImGen package allows you to perform Principal Component Analysis (PCA), partial correlation, gene co-expression, correlated gene expression analysis. With the assistance of a gene enrichment analysis tool (WebGestalt), the significant genes are examined to obtain higher-order biological processes.  Since increasing the number of parcels can inflate statistical significance, NImGen provides a permutation-testing approach based on BrainSMASH. An efficient, reproducible, and rapid pipeline with several MRI markers may be constructed due to the nimgen package's modular design and functionality. The NImGen toolbox includes a set of commonly used tools accompanied by validated, ready-to-use, and easy-to-follow code snippets. These code snippets will guide you through each analysis.

## Requirements

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


## FileTree
```
├── allen_data                                  : gene expression data from Allen Human Brain
├── nimgen                                      : package
├── code
│   ├── compute_gene_lists
│   │   ├── create_submit_files.py              : creates submit files for  ht condor with given parameters
│   │   ├── job1.py                             : distance matrix creation for BrainSmash
│   │   ├── job2.py                             : create surrogate map and  performs correlation analysis (smap x marker)
│   │   ├── job3.py                             : GSE analysis (webgestalt),  coexpression,  region correlation
│   │   ├── job4.py                             : performs PCA
│   │   ├── run_in_venv.sh
│   │   └── test.py
├── input
│   ├── markers
│   │   ├── AOMIC
│   │   └── HCP
│   │       ├──  MarkerDirectory
│   └── parcellations
│       ├── Schaefer2018_100P...                : atlas file directory
│           ├── smaps                           : surrogate map directory (contains surrogate maps for given atlas ie.1000)
│           ├── brain_map.txt                   : needed for BrainSmash   
│           ├── voxel_coordinates.txt           : needed for BrainSmash
│           ├── index.npy                       : needed for BrainSmash
│           ├── distmat.npy                     : needed for BrainSmash
│           ├── Schaefer2018_100P....nii.gz     : parcelleation file       
│           └── nimgen_expressions.csv          : gene expression values for given parcellation. (region X gene)
├── output
│   └── HCP
│       └── Schaefer2018_100Pa..
│           └── VBM
│               └── S1200_AverageT1wDividedByT2w
│                   ├──  smap_corr_scores       : correlation scores for each surrogate map and marker (surrogate map X marker)
│                   ├──  pca_comps              : pca components for original atlas
│                   ├──  Project_genes          : Webgestalt Result
│                   ├── gene_by_gene_matrix.csv : gene by gene correlation score for significants
│                   ├── genes.txt               : significant genes after empirical p-val calculation
│                   ├── pca_comps               : pca components for original atlas
│                   └── reg_by_reg_matrix.csv   : region by region correlation score for significants
└── submit_files                                : submit files for htcondor
```


## Installing

NImGen was developed and tested using Python 3.9. Backwards-compatibility with other Python3 versions is expected but not guaranteed.
NImGen is most easily installed via GitHub repository


```
git clone https://github.com/juaml/nimgen.git
cd nimgen
pip install -e .
```

## Usage

Two components are needed to use the NImGen functionality:

`Atlas/Parcellation`: Brain atlases i.e. Schaefer atlas 100, 500 or 1000. Supported file type: .nii (nifti file).    
`MRI Marker`: fMRI or dMRI marker/weights. Nifti of voxel based morphometry as e.g. outputted by CAT. Extracts measures of region-wise gray matter density (GMD) with mean aggregation method. Supported file type: .nii (nifti file).  
  
Once these are ready you can create surrogate maps for given parcellation file, fetch gene expression data for surrogate maps, calculate empirical p-value and find significant genes.

If you want to use pipeline at once, you can use the command below.  

```
./run_analysis_allinone.py PARCEL_FILE MARKER_FILE N_SURROGATE_MAPS=1000 PARTIAL_CORR=True PERFORM_PCA=True N_PCA_COMP=5
```

However, it is recommended to use seperated jobs. Jobs are separated to four single job to use of computer resources in a most efficient way.

`job1.py`                             : distance matrix creation for BrainSmash  
`job2.py`                             : create single surrogate map for given ID (i.e.: 0_smap.nii) and  performs correlation analysis. This job should be repeated for each surrogate map.  
`job3.py`                            : exports significant genes, GSE analysis (webgestalt), gene coexpression analysis,  region correlation analysis  
`job4.py`                             : performs PCA  

```
./job1.py PARCEL_FILE MARKER_FILE N_SURROGATE_MAPS=1000
./job2.py PARCEL_FILE MARKER_FILE 0
./job3.py PARCEL_FILE MARKER_FILE 
./job4.py PARCEL_FILE MARKER_FILE partial_correlation=True perform_pca=True n_comps=5 
```

for detailed information please check `examples` folder.  

## Creation of submit files for `HT_CONDOR`    

Edit `parcel_and_markers variable` in `create_submit_file.py` and run file.  
It will create 4 seperated condor_submit file.

```
parcels_and_markers = [
        ("input/parcellations/PARCEL_DIR/PARCEL_NAME.nii", "input/markers/HCP/FC/MARKER.nii"),
        ("input/parcellations/PARCEL_DIR/PARCEL_NAME.nii",    "input/markers/HCP/VBM/MARKER.nii"),
    ]
```

## Output Structure

###### Parcellation specific (project_dir/input/parcellations/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm/)  
`distmat.npy & index.npy`: distance-matrix-related memory-mapped files for volumetric data. These files are created by `BrainSmash` package with  `brainsmash.workbench.geo.volume` function and used as inputs to `brainsmash.mapgen.sampled.Sampled` function. File size is ~140 GB.   
`brain_map.txt`: Indicates parcel number of each voxel. Created by export_voxel_coordinates function.  
`voxel_coordinates.txt`: XYZ voxel coordinates from every voxel within an ROI image consisting of ones and zeros based on the parcellation file. Created by export_voxel_coordinates function.  
`smaps`: Contains surrogate maps for given parcellation file and expression file specific for parcellation file. For each surrogate map there should be a unique ID i.e. 1_smap.nii.   

###### Marker specific
`smap_corr_scores`: Contains correlation score files for each surrogate and marker (surrogate map X marker) Created by get_corr_scores function. (function: generate_surrogate_map and get_corr_scores)  
`pca_comps`              : pca components for original atlas (function: get_corr_scores)  
`project_genes`          : Gene set enrichment analysis result, Created by run_webgestalt function. (function: run_webgestalt)  
`genes.txt` : Searches .csv (correlation scores) files in output directory, concats all gene correlation scores of surrogate maps. Exports significant genes based on the all correlation scores for each surrogate map. Applies FDR and export signigicant gene list to genes.txt file. (function: export_significant_genes)  
`pca_comps`               : pca components for original atlas. (function: get_corr_scores)  
`reg_by_reg_matrix.csv`   : region by region correlation score for significant genes. (function: correlated_gene_expression)  
 `gene_by_gene_matrix.csv` : gene by gene correlation score for significant genes. (function: gene_coexpression)  

## Allen Brain Atlas 

Please get the `normalized microarray data` for all the subjects from here:  https://human.brain-map.org/static/download  

Put all of them in a single directory and decompress. The directory should be in the same directory with NImGen package. (check FileTree)  

You will have to provide this directory to the functions calls that get the gene expression data.  

## Core development team
Kaustubh Patil, Fede Raimondo, Leonard Sasse, Talip Yasir Demirtaş  

## Support
If you run into a problem or identify a bug in the source code, please send the authors an email or create a new issue on GitHub.

