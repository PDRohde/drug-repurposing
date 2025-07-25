#------------------------------------------------------------------------------#
# This script performs multi-trait gene-set enrichment analysis using the MAGMA 
# framework. Gene sets based on ATC level 4 drug classifications are retrieved 
# from the gact database. Gene-level Z-scores from VEGAS analyses across nine 
# GWAS traits are combined and analyzed using a Bayesian multi-trait BLR model 
# to identify enriched drug-related pathways.
#------------------------------------------------------------------------------#
# Load libraries
  library(qgg)
  library(gact)

# Load Glist with information on 1000G matched to the ancestry of GWAS data
  Glist <- readRDS(file.path(GAlist$dirs["marker"],"Glist_1000G_eur_filtered.rds"))

# Load GAlist with information on gact database
  GAlist <- readRDS(file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Get gene sets for feature from gact database
  #sets <- getSetsDB(GAlist = GAlist, feature = "DrugGenes", minsets=5)
  #sets <- getSetsDB(GAlist = GAlist, feature = "DrugATCGenes") # Drugs with ATC code
  #sets <- getSetsDB(GAlist = GAlist, feature = "ATC1Genes") # ATC level 1
  #sets <- getSetsDB(GAlist = GAlist, feature = "ATC2Genes") # ATC level 2
  #sets <- getSetsDB(GAlist = GAlist, feature = "ATC3Genes") # ATC level 3
  
  sets <- getSetsDB(GAlist = GAlist, feature = "ATC4Genes") # ATC  level 4

 # Obtain Z-statistics from VEGAS 
  z <- getVEGAS(GAlist,studyID= paste0("GWAS",1:9))
  
# Run Multi-trait BLR model  
  fitMT <- magma(stat=z, sets=sets, method="bayesC", pi=0.01, nit=10000, nburn=1000)
