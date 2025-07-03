#------------------------------------------------------------------------------#
# This script installs the 'gact' and 'qgg' R packages, downloads a versioned 
# gact database of GWAS summary statistics and annotations, and sets up the 
# infrastructure for later use. It also provides commands to explore the 
# structure and content of the downloaded database.
#------------------------------------------------------------------------------#
library(devtools)
devtools::install_github("psoerensen/gact")
devtools::install_github("psoerensen/qgg")

#------------------------------------------------------------------------------#
# Download and install gact database
#------------------------------------------------------------------------------#
# Load libraries
  library(qgg)
  library(gact)

# Define working directory for storing the data base
  dbdir <- "/../projects/gact"

# Create infrastructure and download database (can take several minutes)
  GAlist <- gact(version="hsa.0.0.1", dbdir=dbdir, task="download")

# Save GAlist for use later
  saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Update GAlist for use later
  GAlist <- readRDS(file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Add ATC information to database
  GAlist <- downloadDB(GAlist=GAlist, what="drugbank")
  GAlist <- downloadDB(GAlist=GAlist, what="atc")

  saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

#------------------------------------------------------------------------------#
# Overview of database content
#------------------------------------------------------------------------------#

# Overview of GWAS summary statistics in database
  GAlist$studies

# Accessing the 'dirs' slot in the GAlist object
  GAlist$dirs

# Listing all directories recursively under each path in GAlist$dirs
  list.dirs(GAlist$dirs, recursive = TRUE)

# Listing files under the 'marker' directory specified in GAlist$dirs
  list.files(GAlist$dirs["marker"])

# Listing files under the 'gstat' directory specified in GAlist$dirs
  list.files(GAlist$dirs["gstat"])

# Listing files under the 'gsea' directory specified in GAlist$dirs
  list.files(GAlist$dirs["gsea"])



