##############################################################################
###                                                                        ###
###              Differential expression testing using DESeq2              ###
###                                                                        ###
##############################################################################
#deseq2.r


# remind R where to look for libraries
.libPaths(c("C:/Users/ale097/Data School/Packages"))
# load libraries
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")

rm(list=ls())

df <- read.csv("data/ilrun_counts_siRNA.csv", header=TRUE, row.names=1)
md <- read.csv("data/ilrun_metadata_siRNA.csv", header=TRUE, row.names=1)

all(rownames(md) == colnames (df))

dim(df)

# Calculate counts per million.
# Filter rows: at least 3 samples must have at least 1 cpm.
# Retain rows according to 'keep' and nothing on columns.
cpm <- apply(df, 2, function(x) (x/sum(x))*1000000)
keep <- rowSums( cpm >= 1 ) >=6
df <- df[ keep, ]

# Construct a SummarizedExperiment object:
dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = md,
  design = ~ Condition + Lane) # ~ is representative of 'by', i.e. compare by condition, sex.

# Perform DE testing:
dds <- DESeq(dds)

# Output normalized counts:
norm_counts <- counts (dds, normalized=TRUE)
write.csv(norm_counts, file="Results/DESeq2_ilrun_siRNA_all_normalized_counts.csv")

# Convert results to dataframe:
Inf_NEG24vInf_ILRUN24 <- results(dds, contrast=c("Condition", "Infected_NEG_24hr", "Infected_ILRUN_24hr"))
Uninf_NEG24vUninf_ILRUN24 <- results( dds, contrast=c("Condition", "Uninfected_NEG_24hr", "Uninfected_ILRUN_24hr"))
Uninf_NEG6vInf_NEG6 <- results(dds, contrast=c("Condition", "Uninfected_NEG_6hr", "Infected_NEG_6hr"))
Uninf_NEG24vInf_NEG24 <- results(dds, contrast=c("Condition", "Uninfected_NEG_24hr", "Infected_NEG_24hr"))
Uninf_ILRUN6vInf_ILRUN6 <- results(dds, contrast=c("Condition", "Uninfected_ILRUN_6hr", "Infected_ILRUN_6hr"))
Uninf_ILRUN24vInf_ILRUN24 <- results(dds, contrast=c("Condition", "Uninfected_ILRUN_24hr", "Infected_ILRUN_24hr"))

# Set adjusted p-value significance (padj) threshold:
alpha <- c( 0.05 )
# Set log2FoldChange threshold:
# As it's log2, > 1 is actually equal to 2-fold change or above.
beta <- c( 1 )
# Set baseMean threshold:
gamma <- c( 10 )

# Filter Infectected siNEG versus siILRUN:
# 'Which' provides positions of 'TRUE' values.
sigInf_NEG24vInf_ILRUN24 <- Inf_NEG24vInf_ILRUN24[ which( Inf_NEG24vInf_ILRUN24$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigInf_NEG24vInf_ILRUN24 <- sigInf_NEG24vInf_ILRUN24[ which( abs(sigInf_NEG24vInf_ILRUN24$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigInf_NEG24vInf_ILRUN24 <- sigInf_NEG24vInf_ILRUN24[ which(sigInf_NEG24vInf_ILRUN24$baseMean > gamma), ]
write.csv(sigInf_NEG24vInf_ILRUN24, file="results/sigInf_NEG24vInf_ILRUN24.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
Inf_NEG24vInf_ILRUN24_hits <- rownames(sigInf_NEG24vInf_ILRUN24)
length(Inf_NEG24vInf_ILRUN24_hits)
write.csv(norm_counts[Inf_NEG24vInf_ILRUN24_hits, ], file="results/DESeq2_sig_Inf_NEG24vInf_ILRUN24_normalized_counts.csv")

# Filter Uninfectected siNEG versus siILRUN:
# 'Which' provides positions of 'TRUE' values.
sigUninf_NEG24vUninf_ILRUN24 <- Uninf_NEG24vUninf_ILRUN24[ which( Uninf_NEG24vUninf_ILRUN24$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigUninf_NEG24vUninf_ILRUN24 <- sigUninf_NEG24vUninf_ILRUN24[ which( abs(sigUninf_NEG24vUninf_ILRUN24$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigUninf_NEG24vUninf_ILRUN24 <- sigUninf_NEG24vUninf_ILRUN24[ which(sigUninf_NEG24vUninf_ILRUN24$baseMean > gamma), ]
write.csv(sigUninf_NEG24vUninf_ILRUN24, file="results/sigUninf_NEG24vUninf_ILRUN24.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
Uninf_NEG24vUninf_ILRUN24_hits <- rownames(sigUninf_NEG24vUninf_ILRUN24)
length(Uninf_NEG24vUninf_ILRUN24_hits)
write.csv(norm_counts[Uninf_NEG24vUninf_ILRUN24_hits, ], file="results/DESeq2_sig_Uninf_NEG24vUninf_ILRUN24_normalized_counts.csv")


# Filter Uninfectected siNEG versus Infected siNEG 6hr:
# 'Which' provides positions of 'TRUE' values.
sigUninf_NEG6vInf_NEG6 <- Uninf_NEG6vInf_NEG6[ which( Uninf_NEG6vInf_NEG6$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigUninf_NEG6vInf_NEG6 <- sigUninf_NEG6vInf_NEG6[ which( abs(sigUninf_NEG6vInf_NEG6$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigUninf_NEG6vInf_NEG6 <- sigUninf_NEG6vInf_NEG6[ which(sigUninf_NEG6vInf_NEG6$baseMean > gamma), ]
write.csv(sigUninf_NEG6vInf_NEG6, file="results/sigUninf_NEG6vInf_NEG6.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
Uninf_NEG6vInf_NEG6_hits <- rownames(sigUninf_NEG6vInf_NEG6)
length(Uninf_NEG6vInf_NEG6_hits)
write.csv(norm_counts[Uninf_NEG6vInf_NEG6_hits, ], file="results/DESeq2_sig_Uninf_NEG6vInf_NEG6_normalized_counts.csv")


# Filter Uninfectected siNEG versus Infected siNEG 24hr:
# 'Which' provides positions of 'TRUE' values.
sigUninf_NEG24vInf_NEG24 <- Uninf_NEG24vInf_NEG24[ which( Uninf_NEG24vInf_NEG24$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigUninf_NEG24vInf_NEG24 <- sigUninf_NEG24vInf_NEG24[ which( abs(sigUninf_NEG24vInf_NEG24$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigUninf_NEG24vInf_NEG24 <- sigUninf_NEG24vInf_NEG24[ which(sigUninf_NEG24vInf_NEG24$baseMean > gamma), ]
write.csv(sigUninf_NEG24vInf_NEG24, file="results/sigUninf_NEG24vInf_NEG24.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
Uninf_NEG24vInf_NEG24_hits <- rownames(sigUninf_NEG24vInf_NEG24)
length(Uninf_NEG24vInf_NEG24_hits)
write.csv(norm_counts[Uninf_NEG24vInf_NEG24_hits, ], file="results/DESeq2_sig_Uninf_NEG24vInf_NEG24_normalized_counts.csv")

# Filter Uninfectected siILRUN versus Infected siILRUN 6hr:
# 'Which' provides positions of 'TRUE' values.
sigUninf_ILRUN6vInf_ILRUN6 <- Uninf_ILRUN6vInf_ILRUN6[ which( Uninf_ILRUN6vInf_ILRUN6$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigUninf_ILRUN6vInf_ILRUN6 <- sigUninf_ILRUN6vInf_ILRUN6[ which( abs(sigUninf_ILRUN6vInf_ILRUN6$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigUninf_ILRUN6vInf_ILRUN6 <- sigUninf_ILRUN6vInf_ILRUN6[ which(sigUninf_ILRUN6vInf_ILRUN6$baseMean > gamma), ]
write.csv(sigUninf_ILRUN6vInf_ILRUN6, file="results/sigUninf_ILRUN6vInf_ILRUN6.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
Uninf_ILRUN6vInf_ILRUN6_hits <- rownames(sigUninf_ILRUN6vInf_ILRUN6)
length(Uninf_ILRUN6vInf_ILRUN6_hits)
write.csv(norm_counts[Uninf_ILRUN6vInf_ILRUN6_hits, ], file="results/DESeq2_sig_Uninf_ILRUN6vInf_ILRUN6_normalized_counts.csv")


# Filter Uninfectected siILRUN versus Infected siILRUN 24hr:
# 'Which' provides positions of 'TRUE' values.
sigUninf_ILRUN24vInf_ILRUN24 <- Uninf_ILRUN24vInf_ILRUN24[ which( Uninf_ILRUN24vInf_ILRUN24$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigUninf_ILRUN24vInf_ILRUN24 <- sigUninf_ILRUN24vInf_ILRUN24[ which( abs(sigUninf_ILRUN24vInf_ILRUN24$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigUninf_ILRUN24vInf_ILRUN24 <- sigUninf_ILRUN24vInf_ILRUN24[ which(sigUninf_ILRUN24vInf_ILRUN24$baseMean > gamma), ]
write.csv(sigUninf_ILRUN24vInf_ILRUN24, file="results/sigUninf_ILRUN24vInf_ILRUN24.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
Uninf_ILRUN24vInf_ILRUN24_hits <- rownames(sigUninf_ILRUN24vInf_ILRUN24)
length(Uninf_ILRUN24vInf_ILRUN24_hits)
write.csv(norm_counts[Uninf_ILRUN24vInf_ILRUN24_hits, ], file="results/DESeq2_sig_Uninf_ILRUN24vInf_ILRUN24_normalized_counts.csv")



#PCA plot for all normalized counts

library(tidyverse)

#read in the metadata
all_metadata <- read_csv("data/ilrun_metadata_siRNA.csv") %>% 
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
all_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts.csv") %>%
  rename(gene_locus = X1) %>% 
  gather(Sample, expression, -gene_locus) %>%
  left_join(all_metadata, by = "Sample")

#scale gene expression for all samples 
scaled_genes <- all_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
pca_genes <- prcomp(scaled_genes)

#tidy data frame for plotting
PCA_data <- pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(all_metadata, by = "Sample") %>%
  spread(PC, expression)

#plot to examine variance for all samples 
text_plot_Name <- ggplot(PCA_data, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = Name), size = 4)+
  labs(title = "LT07 is an outlier")
  
ggsave(filename = "results/PCA_all_samples.png", plot = text_plot_Name, width = 12, height = 10, dpi = 300, units = "cm")
#LT07 is an outlier

#####Repeat analysis for dataset excluding LT07###############
#read in the metadata
lessLT07_metadata <- read_csv("data/ilrun_metadata_siRNA.csv") %>% 
  mutate(Name = str_extract(Sample, "^....")) %>%
  filter(Sample != "LT07_L001" |Sample != "LT07_L002")
  

#join with the counts data 
lessLT07_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts.csv") %>%
  rename(gene_locus = X1) %>%
  select(-LT07_L001, -LT07_L002) %>% 
  gather(Sample, expression, -gene_locus) %>%
  left_join(lessLT07_metadata, by = "Sample")

#scale gene expression for all samples 
lessLT07_scaled_genes <- lessLT07_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
lessLT07_pca_genes <- prcomp(lessLT07_scaled_genes)

#tidy data frame for plotting
PCA_data_lessLT07 <- lessLT07_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(all_metadata, by = "Sample") %>%
  spread(PC, expression)


text_plot_lessLT07 <- ggplot(PCA_data_lessLT07, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = siRNA, color = Infection), size = 3)+
  labs(title = "Principle Component Analysis less LT07")
  
ggsave(filename = "results/PCA_lessLT07.png", plot = text_plot_lessLT07, width = 12, height = 10, dpi = 300, units = "cm")
