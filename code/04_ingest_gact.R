#------------------------------------------------------------------------------#
# This script demonstrates how to ingest multiple GWAS summary statistics datasets 
# into the gact database. For each dataset, quality control is performed and 
# allele information is matched against the gact reference. The formatted data 
# is then added to the database using `updateStatDB()`, along with metadata 
# describing trait type, ancestry, sample size, and publication reference.
#------------------------------------------------------------------------------#
# Load libraries
  library(qgg)
  library(gact)
  library(data.table)

# Load GAlist
  GAlist <- readRDS(file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

#------------------------------------------------------------------------------#
# Cardiogram CAD GWAS
#------------------------------------------------------------------------------#
# link: https://cardiogramplusc4d.org/data-downloads/


#------------------------------------------------------------------------------#
# T2D 
#------------------------------------------------------------------------------#
# link: https://www.diagram-consortium.org/downloads.html



#------------------------------------------------------------------------------#
# UK Biobank 2021 HbA1c GWAS: European ancestry
#------------------------------------------------------------------------------#
# link: https://t2d.hugeamp.org/dinspector.html?dataset=SinnottArmstrong2021_HBA1C_EU

# Load GWAS data
  fname_stat <- "/../projects/gact/hsa.0.0.1/download/Glycated_haemoglobin_HbA1c.imp.gz"
  stat <- fread(fname_stat, data.table = FALSE)
  
# Subset and rename columns according to required format
  stat <- stat[, c("MarkerName", "#CHROM", "POS", "ALT", "REF", "Effect", "StdErr", "P-value")]
  colnames(stat) <- c("marker", "chr", "pos", "ea", "nea", "b", "seb", "p")

# Update database
  GAlist <- updateStatDB(GAlist = GAlist,
                         stat = stat,
                         source = "Glycated_haemoglobin_HbA1c.imp",
                         trait = "HbA1C",
                         type = "quantitative",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID:33462484",
                         n = 318779,
                         ncase = 0,
                         ncontrol = 0,
                         comments = "Just UK biobank",
                         writeStatDB = TRUE)

# Save updated database
  saveRDS(GAlist, file = "../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

#------------------------------------------------------------------------------#
# CKDGen 2019 GWAS
#------------------------------------------------------------------------------#
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_CKDGenConsortium

# Load GWAS data
  fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"
  stat <- fread(fname_stat, data.table=FALSE)

# Subset and rename columns according to required format
  stat <- stat[, c("RSID","Chr","Pos_b37","Allele1","Allele2","Freq1","Effect", "StdErr", "P-value", "n_total_sum")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
  GAlist <- updateStatDB(GAlist=GAlist,
                         stat=stat,
                         source="CKD_overall_EA_JW_20180223_nstud23.dbgap.txt",
                         trait="CKD",
                         type = "binary",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID: 31152163",
                         n = 480698,
                         ncase = 41395,
                         ncontrol = 439303,
                         comments ="Include UK biobank",
                         writeStatDB=TRUE)

saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

#------------------------------------------------------------------------------#
# Hypertension GWAS
#------------------------------------------------------------------------------#
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Zhu2019_COPD_CVD_eu

# Load GWAS data
  fname_stat <- "/../projects/gact/hsa.0.0.1/download/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.doctor_highbloodpressure.assoc.gz"
  stat <- fread(fname_stat, data.table=FALSE)

# Subset and rename columns according to required format
  stat <- stat[, c("SNP","CHR","BP","A1","A0", "BETA", "SE", "P", "INFO")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "b", "seb", "p", "info")

# Update database
  GAlist <- updateStatDB(GAlist=GAlist,
                         stat=stat,
                         source="ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.doctor_highbloodpressure.assoc",
                         trait="HTN",
                         type = "binary",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID:30940143",
                         n = 458554,
                         ncase = 144793,
                         ncontrol = 313761,
                         comments ="Include UK biobank",
                         writeStatDB=TRUE)

saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

#------------------------------------------------------------------------------#
#  Meta-analysis of body mass index (bmi) in UK Biobank and GIANT data.
#------------------------------------------------------------------------------#
# link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6298238/

# Load GWAS data
  fname_stat <- "/../projects/gact/hsa.0.0.1/download/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
  stat <- fread(fname_stat, data.table=FALSE)

# Subset and rename columns according to required format
  stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele", "Freq_Tested_Allele","BETA", "SE", "P", "N","INFO")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "info")

# Update database
  GAlist <- updateStatDB(GAlist=GAlist,
                         stat=stat,
                         source="Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt",
                         trait="BMI",
                         type = "quantitative",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID:30239722",
                         n = 806834,
                         ncase = 0,
                         ncontrol = 806834,
                         comments ="Include UK biobank",
                         writeStatDB=TRUE)
  saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

#------------------------------------------------------------------------------#
# Meta-analysis of whr in UK Biobank and GIANT data.
#------------------------------------------------------------------------------#
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_UKBiobankGIANT_eu

# Load GWAS data
  fname_stat <- "/../projects/gact/hsa.0.0.1/download/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
  stat <- fread(fname_stat, data.table=FALSE)

# Subset and rename columns according to required format
  stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele", "Freq_Tested_Allele","BETA", "SE", "P", "N","INFO")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "info")

# Update database
  GAlist <- updateStatDB(GAlist=GAlist,
                         stat=stat,
                         source="Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt",
                         trait="WHR",
                         type = "quantitative",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID:30239722",
                         n = 694649,
                         ncase = 0,
                         ncontrol = 694649,
                         comments ="Include UK biobank",
                         writeStatDB=TRUE)

  saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

#------------------------------------------------------------------------------#
# Triglycerides
#------------------------------------------------------------------------------#
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Graham_2021_lipids_Mixed

# Load GWAS data
  fname_stat <- "/../projects/gact/hsa.0.0.1/download/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
  stat <- fread(fname_stat, data.table=FALSE)

# Subset and rename columns according to required format
  stat <- stat[, c("rsID","CHROM","POS_b37","ALT","REF", "POOLED_ALT_AF", "EFFECT_SIZE", "SE", "pvalue", "N")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
  GAlist <- updateStatDB(GAlist=GAlist,
                         stat=stat,
                         source="logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
                         trait="TG",
                         type = "quantitative",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID:34887591",
                         n = 1320016,
                         ncase = 0,
                         ncontrol = 1320016,
                         comments ="Include UK biobank",
                         writeStatDB=TRUE)

  saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

#------------------------------------------------------------------------------#
# Systolic blood presure
#------------------------------------------------------------------------------#
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Evangelou2018_bp_eu

# Load GWAS data
  fname_stat <- "/../projects/gact/hsa.0.0.1/download/Evangelou_30224653_SBP.txt.gz"
  stat <- fread(fname_stat, data.table=FALSE)

# Modify columns according to required format
  cp <- strsplit(stat$MarkerName, ":")
  stat$chr <- sapply(cp, function(x){x[1]})
  stat$pos <- as.numeric(sapply(cp, function(x){x[2]}))
  stat$Allele1 <- toupper(stat$Allele1)
  stat$Allele2 <- toupper(stat$Allele2)

# Subset and rename columns according to required format
  stat <- stat[, c("MarkerName", "chr","pos", "Allele1","Allele2","Freq1","Effect", "StdErr", "P", "N_effective")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
  GAlist <- updateStatDB(GAlist=GAlist,
                         stat=stat,
                         source="Evangelou_30224653_SBP.txt",
                         trait="SBP",
                         type = "quantitative",
                         gender = "both",
                         ancestry = "EUR",
                         build = "GRCh37",
                         reference = "PMID:302246538",
                         n = 757601,
                         ncase = 0,
                         ncontrol = 757601,
                         comments ="Include UK biobank",
                         writeStatDB=TRUE)

  saveRDS(GAlist, file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)