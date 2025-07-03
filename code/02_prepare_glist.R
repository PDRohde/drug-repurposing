#------------------------------------------------------------------------------#
# This script prepares genotype data summaries (Glist) for the 1000 Genomes 
# Project (1000G) within the gact framework. It loads required libraries, 
# downloads the 1000G reference data if needed, and defines paths for genotype 
# files. The script creates a genotype summary object for the entire 1000G 
# dataset, and also processes European (EUR) ancestry data specifically by 
# filtering variants based on those available in the gact database. 
# Filtered genotype files are saved and prepared as separate Glist objects for 
# downstream analyses. Similar procedures can be applied for other ancestries as well.
#------------------------------------------------------------------------------#

# Load libraries
  library(qgg)
  library(gact)
  library(data.table)

# Load GAlist
  GAlist <- readRDS(file="../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Download 1000G data (if not allready downloaded)
  GAlist <- downloadDB(GAlist=GAlist, what="1000G")

# Define location of bed/bim/fam files
  path <- file.path(GAlist$dirs["marker"],"1000G_EUR_Phase3_plink/1000G.EUR.QC.")
  bedfiles <- paste(path,1:22,".bed",sep="")
  bimfiles <- paste(path,1:22,".bim",sep="")
  famfiles <- paste(path,1:22,".fam",sep="")

# Prepare summary (i.e. Glist) for 1000G genotype data
  Glist <- gprep(study="1000G",
                 bedfiles=bedfiles,
                 bimfiles=bimfiles,
                 famfiles=famfiles)

# Save Glist for use later
  saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G.rds"))


#------------------------------------------------------------------------------#
# Prepare Glist for 1000G data for different ancestries
#------------------------------------------------------------------------------#

# Load GAlist
   GAlist <- readRDS(file="../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Marker IDs in database
  rsids <- GAlist$rsids


#------------------------------------------------------------------------------#
# Process EUR 1000G data
#------------------------------------------------------------------------------#

# Define the file paths for the original bed/bim/fam files to read
  bedRead <- file.path(GAlist$dirs["marker"], "g1000_eur.bed")
  bimRead <- file.path(GAlist$dirs["marker"], "g1000_eur.bim")
  famRead <- file.path(GAlist$dirs["marker"], "g1000_eur.fam")

# Define the file paths for the filtered bed/bim/fam files to write
  bedWrite <- file.path(GAlist$dirs["marker"], "g1000_eur_filtered.bed")
  bimWrite <- file.path(GAlist$dirs["marker"], "g1000_eur_filtered.bim")
  famWrite <- file.path(GAlist$dirs["marker"], "g1000_eur_filtered.fam")

# Call the writeBED function to filter and write the data
  writeBED(bedRead = bedRead,
           bimRead = bimRead,
           famRead = famRead,
           bedWrite = bedWrite,
           bimWrite = bimWrite,
           famWrite = famWrite,
           rsids = rsids
  )

# Define location of bed/bim/fam files
  bedfiles <- file.path(GAlist$dirs["marker"],"g1000_eur_filtered.bed")
  bimfiles <- file.path(GAlist$dirs["marker"],"g1000_eur_filtered.bim")
  famfiles <- file.path(GAlist$dirs["marker"],"g1000_eur_filtered.fam")

# Prepare summary (i.e. Glist) for 1000G genotype data
  Glist <- gprep(study="1000G EUR",
                 bedfiles=bedfiles,
                 bimfiles=bimfiles,
                 famfiles=famfiles)

# Save Glist for use later
  saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G_eur_filtered.rds"))

# Similar can be done for the other ancestries