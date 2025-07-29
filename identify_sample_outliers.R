message('Script Begin')

library(RNOmni)
library(tidyverse)
library(data.table)
library(WGCNA)
library(optparse)
library(R.utils)
####### COMMAND LINE ARGUMENTS ########

message('Parsing options')
option_list <- list(
    optparse::make_option(c("--TPM_file"), type="character", default=NULL,
                        help="Sample to be used in processing expression marix", metavar = "type"),
    optparse::make_option(c("--count_file"), type="character", default=NULL,
                        help="sample to be used in processing expression marix", metavar = "type"),
    optparse::make_option(c("--prefix"), type="character", default=NULL,
                        help="sample to be used in processing expression marix", metavar = "type")
    )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

message('Loading command line arguments')
RNAseqQC2_TPMs <- opt$TPM_file 
RNAseqQC2_counts <- opt$count_file 

outliers_file <- paste0(opt$prefix,'_connectivity_outliers.tsv')
all_z_scores_file <- paste0(opt$prefix,'_connectivity_scores.tsv')


####### BEGIN #########

message('Loading counts and TPMs')
TPM_df <- fread(RNAseqQC2_TPMs,skip = 2,header= TRUE)
counts_df <- fread(RNAseqQC2_counts,skip = 2,header = TRUE)

transposed_counts <- counts_df %>% 
    select(-Description) %>% 
    column_to_rownames('Name') %>% 
    t() %>% 
    data.frame()

expressed_genes <- transposed_counts %>% 
    dplyr::select(where(~ mean(.x > 6) >= 0.2))

message('Normalizing TPMs')
rotated_normalized_TPMs <- TPM_df %>% 
    filter(Name %in% colnames(expressed_genes)) %>% 
    select(-Description) %>% 
    column_to_rownames('Name') %>% 
    t() %>% 
    data.frame() %>% 
    mutate(across(everything(),~RankNorm(.)))


##### COMPUTE OUTLIERS USING WGCNA ######

# compute connectivity score by calculating correlation 
# of all samples with each other and then Z-scoring the data
message('Computing correlation matrix')
norm_adj <- (0.5 + 0.5 * bicor(rotated_normalized_TPMs))

message('Computing network')
net_summary <- fundamentalNetworkConcepts(norm_adj)
net_connectivity <- net_summary$Connectivity
connectivity_zscore <- ((net_connectivity - mean(net_connectivity)) / sd(net_connectivity)) %>% 
        data.frame()  %>% 
        dplyr::rename('Z_score' = 1) %>% 
        rownames_to_column('ResearchID')

connectivity_zscore_outliers <- connectivity_zscore %>% filter(Z_score < -3)

message('Writing data')
connectivity_zscore_outliers %>% write_tsv(outliers_file)
connectivity_zscore %>% write_tsv(all_z_scores_file)


