#------------------------------------------------------------------------------#
# This script computes sparse linkage disequilibrium (LD) matrices for the 1000 
# Genomes Project (1000G) genotype data using the gact framework. It loads the 
# necessary libraries and previously prepared genotype summaries (Glist). 
# The script identifies high-quality genetic markers based on several quality 
# filters (e.g., minor allele frequency, missingness, indels, duplicates) and 
# computes sparse LD matrices for these markers. The resulting Glist is updated
# and saved. Additionally, a marker annotation table containing positions, 
# alleles, allele frequencies, and LD scores is created and exported for 
# downstream analyses. For efficient computation, linking R to optimized BLAS 
# libraries (e.g., OpenBLAS, MKL) is recommended.
#------------------------------------------------------------------------------#

# Load libraries
  library(qgg)
  library(gact)
  library(data.table)

# Load GAlist
  GAlist <- readRDS(file="/../projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Load Glist 
  Glist <- readRDS(file=file.path(GAlist$dirs["marker"],"Glist_1000G.rds"))


#------------------------------------------------------------------------------#
# Compute sparse LD matrices used for computing LD matrices in downstream analyses
#------------------------------------------------------------------------------#
# This step may take some tome to complete. In genereal it is recommended that R
# is linked to openblas, MKL or Atlas (default MacOS) for fast computation

# Identify high quality genetic markers in 1000G used for computation of sparse LD matrices
  rsids <- gfilter(Glist = Glist,
                   excludeMAF = 0.05,
                   excludeMISS = 0.05,
                   excludeCGAT = TRUE,
                   excludeINDEL = TRUE,
                   excludeDUPS = TRUE,
                   excludeHWE = 1e-12,
                   excludeMHC = FALSE)

  Glist <- gprep(Glist, task = "sparseld", msize = 1000, rsids = rsids, overwrite = TRUE)
  saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G.rds"))

  chr <- 1:22
  markers <- data.frame(rsids=unlist(Glist$rsids[chr]),
                        chr=unlist(Glist$chr[chr]),
                        pos=unlist(Glist$pos[chr]),
                        ea=unlist(Glist$a1[chr]),
                        nea=unlist(Glist$a2[chr]),
                        eaf=unlist(Glist$af[chr]),
                        maf=unlist(Glist$maf[chr]),
                        map=unlist(Glist$map[chr]),
                        ldscores=NA)
  rownames(markers) <- markers$rsids
  markers[names(unlist(Glist$ldscores[chr])),"ldscores"] <- unlist(Glist$ldscores[chr])
  fwrite(markers, file=file.path(GAlist$dirs["marker"],"markers_1000G.txt.gz"))

