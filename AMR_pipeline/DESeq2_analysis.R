################################### Install Packages ################################### 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
########################################################################################


######################################## STEP 1 ########################################
# Load the AMR_meta.txt file from AMR_pipeline
amr_meta = read.csv("/path/to/AMR_pipeline/AMR_meta.txt", sep="\t")

# Load your metadata file containing sample IDs and treatment conditions
# Sample IDs must match the sample IDs in your Salmon qunatification output
# There can only be two treatment conditions
conditions = read.csv("/path/to/conditions/file.csv", row.names=1)

# Load Salmon results, saving NumReads in one file
# There should only be the desired output files in salmon_results_dir
salmon_results_dir = "/path/to/salmon/csv/results/" # Directory your csv files are
salmon_files = list.files(salmon_results_dir) # List all the files in the directory

# Initiate an empty dataframe to save all NumReads columns in
all_reads = data.frame(row.names = amr_meta$MEG_ID)

# Iterate through each file, selecting the NumReads and concatenating to all_reads
# Column names are made sample ID
for (file in salmon_files) {
  
  sample = gsub(".tsv", "", file)
  tsv = read.csv(paste0(salmon_results_dir, file), sep="\t") %>% 
    select(NumReads)
  colnames(tsv) = c(sample)
  
  all_reads = cbind(all_reads, tsv[sample])
}
########################################################################################


######################################## STEP 2 ########################################
# Get the name of the shared IDs between conditions and all_reads
samples = intersect(colnames(all_reads), row.names(conditions))

# Select just the shared samples in the meta data
conditions_filter = conditions[row.names(conditions) %in% samples,]

# Ensure samples are in the same order in all_reads
all_reads = all_reads %>% select(row.names(conditions_filter))

# Prep DESeqDataSet object
conditions_filter$column_name = factor(conditions_filter$column_name,         # Change this to the desired column CHANGE 2 LOCATIONS
                                           levels= c("Level 1", "Level 2"))   # Set factors CHANGE

all_reads <- round(all_reads,0)# Round count data
########################################################################################


######################################## STEP 3 ########################################
# Create DESeqDataSet (DO NOT EDIT THIS)
dds <- DESeqDataSetFromMatrix(countData = all_reads, # NumReads dataframe, rounded to nearest integer
                              colData = conditions_filter, # Conditions dataframe
                              design = ~cond_col) # Specified conditions column in variable cond_col CHANGE

# Perform DGE analysis
dds <- DESeq(dds)
########################################################################################


######################################## STEP 4 ########################################
# Get unfiltered results
res <- results(dds) #unfiltered
write.csv(res, "/path/to/save/unfiltered_results.csv")
# head(res) # Use this to view the first 5 rows




# Get filtered results
pvalue = 0.05 # pvalue threshold to filter results, less stringent than padj
log2FC = 1 # The absolute log2FC value must be greater than this

res_filt <- res[which(res$pvalue < pvalue & abs(res$log2FoldChange) > log2FC), ]
write.csv(res_filt, "/path/to/save/filtered_results.csv")
# head(res_filt) # Use this to view the first 5 rows
########################################################################################


##################################### SNP Analysis #####################################
cond_col = "condition_column" # Name of the column containing the desired condition to evaluate
group1 = "condition1"
group2 = "condition2"

# Create a new df with just sample names and condition labels
for_snp_analysis = data.frame(conditions_filter) %>% select(cond_col)
for_snp_analysis$samples = row.names(for_snp_analysis)

# Subset the df into 2 df, one for each group
group1_df = for_snp_analysis[for_snp_analysis[cond_col] == group1,]
group2_df = for_snp_analysis[for_snp_analysis[cond_col] == group2,]

# Change the sample names to the file names to merge together
group1_df$samples = paste0(group1_df$samples, "_variant_filter.vcf.gz")
group2_df$samples = paste0(group2_df$samples, "_variant_filter.vcf.gz")

# Save as lists as text files
write.table(group1_df$samples, "/path/to/save/group1.txt", # File path and name for group 1
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(group2_df$samples, "/path/to/save/group1.txt", # File path and name for group 2
            sep="\t", quote=F, row.names=F, col.names=F)
########################################################################################


