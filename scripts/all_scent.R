all_scent <- function(subgroup, beach) {
# Main script for analysing scent profiles and associations with genetics

# laoding functions ------------------------------------------------------------
source("get_scores.R")
source("subset_all.R")
source("get_pairdiff.R")
source("rel_results.R")
source("het_results.R")
# loading data------------------------------------------------------------------                                                                                                                                                                                                                                                                                                                               

# abundance matrix (standardized and log(x + 1) - transformed)
# alternatives: scent abundances.csv,  scent abundance nobeach.csv,
#               scent.abundance.pa.csv

scent_abundance <- as.data.frame(t(read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                                  "MSc.Behaviour\\Research\\Seal Scent\\R code\\",
                                                  "data\\csv_files\\scent abundances.csv", 
                                                  sep = ""), row.names=1)))

# diversity measures
scent_diversity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "MSc.Behaviour\\Research\\Seal Scent\\",
                                  "R code\\data\\csv_files\\",
                                  "scent diversity.csv", sep = ""),
                            row.names=1)

# relatedness matrix (old: relatedness_41loci.csv)
relatedness <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                              "MSc.Behaviour\\Research\\Seal Scent\\",
                              "R code\\data\\csv_files\\",
                              "relatednessnew.csv", sep = ""),
                        row.names=1)

## heterozygosity SH
heterozygosity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                 "MSc.Behaviour\\Research\\Seal Scent\\",
                                 "R code\\data\\csv_files\\",
                                 "heterozygosity_41loci.csv", sep = ""),
                           row.names=1) 

# beach and family factor
factors <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                          "MSc.Behaviour\\Research\\Seal Scent\\",
                          "R code\\data\\csv_files\\",
                          "factors.csv", sep = ""),
                    row.names=1) 


# choose ordination method and get scores---------------------------------------

allscores <- get_scores(scent_abundance, method = "fa", num_dim = 4,
                        rotation = "promax")

# get subsets for age ("mums", "pups", "all") and beach (0 = all, 1, 2)---------

allsubsets <- subset_all(subgroup, beach, scent_abundance,scent_diversity,
                         relatedness, heterozygosity, allscores, factors)

abund <- allsubsets[[1]]
div <- allsubsets[[2]]
relate <- allsubsets[[3]]
het <- allsubsets[[4]]
scores <- allsubsets[[5]]
factors <- allsubsets[[6]]

num_fac <- ncol(scores)

# relatedness-------------------------------------------------------------------

# compute relatedness + scores data frame
relate_df <- get_pairdiff(relate, scores, df = F)

# and list of relatedness-difference matrices
relate_list <- get_pairdiff(relate, scores, df = T)

# relatedness results
relate_results <- rel_results(relate, scores, abund, num_fac)

# heterozygosity----------------------------------------------------------------

# create main data.frame for heterozygosity 
het_df <- data.frame("het"= het[,1],
                     num_comps = div$S,
                     row.names=rownames(het))
het_df <- cbind(het_df, scores[, 1:ncol(scores)])

# heterozygosity results
het_results <- het_results(het_df, num_fac)

out <- list(relate.df=relate_df,
             results.relatedness=relate_results,
             relate.list=relate_list,
             relatedness=relate,
             het.df=het_df,
             results.het=het_results,
             scores=scores,
             factors=factors)
}