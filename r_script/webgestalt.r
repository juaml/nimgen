#! /usr/bin/Rscript
# -------------------------------------
# NImGen:
# Multiarray analysis of Allen Human Brain Atlas
# Colaborative project with Kaustubh, Leo, and Amir
# -------------------------------------

if (!require("WebGestaltR"))
  install.packages("WebGestaltR")
library(WebGestaltR)
library(tools)

# ------------------------------------------
# ----- Set the initial parameters and folder addresses
# ------------------------------------------
args <- commandArgs(TRUE)
input_file <- args[1]
input_file_path <- file_path_as_absolute(input_file)
file_name <- file_path_sans_ext(basename(input_file))
file_ext <- file_ext(input_file)
output_path <- dirname(input_file_path)

# ------------------------------------------
# ----- parameters for analysis type
# ------------------------------------------

enrichDatabase <- c(
  "disease_GLAD4U",
  "pathway_KEGG",
  "phenotype_Human_Phenotype_Ontology",
  "drug_DrugBank",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function",
  "geneontology_Biological_Process"
)
organism <- "hsapiens"
ref_set <- "agilent_sureprint_g3_ge_8x60k"

# ----------------------------------------------
# Main script
# ----------------------------------------------
if (!is.na(input_file)) {
  tryCatch({
      enrichResult <- WebGestaltR(
        enrichMethod = "ORA",
        organism = organism,
        enrichDatabase = enrichDatabase,
        interestGeneType = "genesymbol",
        interestGeneFile=input_file,
        referenceSet = ref_set,
        referenceGeneType = "genesymbol",
        isOutput = TRUE,
        outputDirectory = output_path,
        projectName = file_name
      )
      print("DONE")
    },
    error = function(e) {
      print("Error: WebGestaltR couldn't get data.")
      print(e)
    }
  )
  
} else {
  print("Error: Some arguments are missing: input_file")
}
