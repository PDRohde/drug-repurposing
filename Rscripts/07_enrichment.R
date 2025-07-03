#-------------------------------------------------------------------------------
# This script performs enrichment analysis to identify drug classes (ATC level 4) 
# that are significantly associated with diabetes-related gene sets. 
# Gene sets are tested using hypergeometric tests based on different sources of 
# disease-gene associations from the gact database, including curated knowledge, 
# text mining, and experimental data.
#-------------------------------------------------------------------------------
# Load libraries
  library(qgg)
  library(gact)

# Load Glist with information on 1000G matched to the ancestry of GWAS data
  Glist <- readRDS(file.path(GAlist$dirs["marker"],"Glist_1000G_eur_filtered.rds"))

# Load GAlist with information on gact database
  GAlist <- readRDS(file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Define gene sets
  sets <- getSetsDB(GAlist = GAlist, feature = "ATC4Genes") # ATC  level 4

# Run hypergeometric test
  res <- hgtDB(GAlist = GAlist, sets = sets, feature = "DiseaseGenes", output="p") 
  rws <- c("N03AF" ,"A10AE","D10AX", "A10BX", "C10AB","L01EF","L04AB", "L01XX","L01XG", "M04AA",
           "L01EE", "P02BA","R06AB", "D06AA","L01XD")
  cls <- grep("diabetes", tolower(colnames(res)))
  head(t(res[rws,cls]), 50)

# Enrichment test based on textmining information
  res <- hgtDB(GAlist = GAlist, sets = sets, feature = "DiseaseGenesTMplus", output="p") 
  cls <- grep("diabetes", tolower(colnames(res)))
  head(t(res[rws,cls]), 50)

# Enrichment test based on textmining information
  res <- hgtDB(GAlist = GAlist, sets = sets, feature = "DiseaseGenesKBplus", output="p") 
  cls <- grep("diabetes", tolower(colnames(res)))
  head(t(res[rws,cls]), 50)

# Enrichment test based on experiment information
  res <- hgtDB(GAlist = GAlist, sets = sets, feature = "DiseaseGenesEXPplus", output="p") 
  cls <- grep("diabetes", tolower(colnames(res)))
  head(t(res[rws,cls]), 50)

# Extract the subset of the matrix
  cls <- grep("diabetes", tolower(colnames(res)))
  subset_matrix <- res[rws, cls]
