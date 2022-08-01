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
analysis_type <- args[2]
if (is.na(analysis_type)) {
  analysis_type = FALSE
}
input_file_path <- file_path_as_absolute(input_file)
file_name <- file_path_sans_ext(basename(input_file))
file_ext <- file_ext(input_file)
output_path <- dirname(input_file_path)


# ------------------------------------------
# ----- read different types of files
# ------------------------------------------

tryCatch( { 
  data = read.delim(input_file, header = FALSE) 
  }
, error = function(e) 
  {
  print("Error occured while reading input file..")

  }
)

# ------------------------------------------
# ----- parameters for analysis type
# ------------------------------------------
if (analysis_type == 'GSEA') {
  geneandscore = data.frame(gene = data[, 1] ,
                            score = format(data[, 3], scientific = FALSE))
  geneandscore$score = as.integer(geneandscore$score)
  enrichDatabase <- "pathway_KEGG"
} else{
  data = data[, 1]
  enrichDatabase <- c(
    "disease_GLAD4U",
    "pathway_KEGG",
    "phenotype_Human_Phenotype_Ontology",
    "drug_DrugBank",
    "geneontology_Cellular_Component",
    "geneontology_Molecular_Function",
    "geneontology_Biological_Process"
  )
}

organism <- "hsapiens"
ref_set <- "agilent_sureprint_g3_ge_8x60k"

# ----------------------------------------------
# Main script
# ----------------------------------------------
if (!is.na(input_file)) {
  tryCatch(
    expr = {
      if (analysis_type != 'GSEA') {
        enrichResult <- WebGestaltR(
          enrichMethod = "ORA",
          organism = organism,
          enrichDatabase = enrichDatabase,
          interestGene = data,
          interestGeneType = "genesymbol",
          # interestGeneFile=input_file_path
          referenceSet = ref_set,
          referenceGeneType = "genesymbol",
          isOutput = TRUE,
          outputDirectory = output_path,
          projectName = file_name
        )
      } else{
        enrichResult <- WebGestaltR(
          enrichMethod = "GSEA",
          organism = "hsapiens",
          enrichDatabase = enrichDatabase,
          interestGene = geneandscore,
          interestGeneType = "genesymbol",
          sigMethod = "top",
          topThr = 10,
          minNum = 5,
          outputDirectory = output_path,
          projectName = paste("GSEA",file_name,sep="_")
        )
      }
      
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
