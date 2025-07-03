#------------------------------------------------------------------------------#
# This script performs gene-level association analysis using the VEGAS method 
# for all GWAS studies available in the gact database. For each study, summary 
# statistics are retrieved, aligned to the 1000G reference panel, and tested 
# for enrichment using predefined gene-marker sets (±40kb/10kb). Results are 
# saved separately for each study.
#------------------------------------------------------------------------------#

# Load libraries
  library(qgg)
  library(gact)

# Load Glist with information on 1000G matched to the ancestry of GWAS data
  Glist <- readRDS(file.path(GAlist$dirs["marker"],"Glist_1000G_eur_filtered.rds"))

# Load GAlist with information on gact database
  GAlist <- readRDS(file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Check studies in gact database
  GAlist$studies

# Extract gene-marker sets (include markers 40kb/10kb upstream/downstream)
  markerSets <- getMarkerSets(GAlist = GAlist, feature = "Genesplus")

# Perform VEGAS gene analysis for a single study of European ancestry (EUR)
for(i in 1:lenght(GAlist$studies)){
	studyID <- paste0("GWAS",i)

	# Get GWAS summary statistics from gact database
	stat <- getMarkerStat(GAlist=GAlist, studyID=studyID)

	# Check and align summary statistics based on marker information in Glist
	stat <- checkStat(Glist=Glist, stat=stat)

	# Gene analysis using VEGAS
	res <- vegas(Glist=Glist, sets=markerSets, stat=stat, verbose=TRUE)
	filename <- file.path(GAlist$dirs["gsea"], paste0(studyID, "_vegas.rds"))
	saveRDS(res,file=filename)
}


