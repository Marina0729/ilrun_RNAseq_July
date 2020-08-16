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
library(DESeq2)
library(tidyverse)
library(dplyr)
library(readr)
library(vroom)
library(stringr)

rm(list=ls())

df <- read.csv("data/ilrun_counts_siRNA.csv", header=TRUE, row.names=1)
md <- read.csv("data/ilrun_metadata_siRNA.csv", header=TRUE, row.names=1)

df_LT07 <- read.csv("data/ilrun_counts_siRNA.csv", header=TRUE, row.names=1) %>% 
  select(-LT07_L001, -LT07_L002)
md_LT07 <- read.csv("data/ilrun_metadata_siRNA_LT07.csv", header=TRUE, row.names=1)

all(rownames(md) == colnames (df))
all(rownames(md_LT07) == colnames (df_LT07))

dim(df)
dim(df_LT07)


# Calculate counts per million.
# Filter rows: at least 3 samples must have at least 1 cpm.
# Retain rows according to 'keep' and nothing on columns.
cpm <- apply(df, 2, function(x) (x/sum(x))*1000000)
keep <- rowSums( cpm >= 1 ) >=4
df_filtered <- df[ keep, ]

cpm_LT07 <- apply(df_LT07, 2, function(x) (x/sum(x))*1000000)
keep_LT07 <- rowSums( cpm_LT07 >= 1 ) >=4
df_filtered_LT07 <- df_LT07[ keep_LT07, ]

dim(df_filtered)
dim(df_filtered_LT07)


# Construct a SummarizedExperiment object:
dds <- DESeqDataSetFromMatrix(
  countData = df_filtered,
  colData = md,
  design = ~ Condition + Lane) # ~ is representative of 'by', i.e. compare by condition

dds_LT07 <- DESeqDataSetFromMatrix(
  countData = df_filtered_LT07,
  colData = md_LT07,
  design = ~ Condition + Lane) # ~ is representative of 'by', i.e. compare by condition

# Perform DE testing:
dds <- DESeq(dds)
dds_LT07 <- DESeq(dds_LT07)

# Output normalized counts:
norm_counts <- counts (dds, normalized=TRUE)
write.csv(norm_counts, file="results/DESeq2_ilrun_siRNA_all_normalized_counts.csv")

norm_counts_LT07 <- counts (dds_LT07, normalized=TRUE)
write.csv(norm_counts_LT07, file="results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07.csv")


# Convert results to dataframe:
CoV2genesNEG6hr <- results(dds, contrast = c("Condition", "Uninfected_NEG_6hr", "Infected_NEG_6hr"))
CoV2genesNEG24hr <- results(dds, contrast = c("Condition", "Uninfected_NEG_24hr", "Infected_NEG_24hr"))
CoV2genesNEG6hr <- results(dds, contrast = c("Condition", "Uninfected_ILRUN_6hr", "Infected_ILRUN_6hr"))
CoV2genesNEG24hr <- results(dds, contrast = c("Condition", "Uninfected_ILRUN_24hr", "Infected_ILRUN_24hr"))
ILRUNgenesUninf6hr <- results(dds, contrast = c("Condition","Uninfected_NEG_6hr", "Uninfected_ILRUN_6hr"))
ILRUNgenesUninf24hr <- results(dds, contrast = c("Condition","Uninfected_NEG_24hr", "Uninfected_ILRUN_24hr"))
ILRUNgenesInf6hr <- results(dds, contrast = c("Condition","Infected_NEG_6hr", "Infected_ILRUN_6hr"))
ILRUNgenesInf24hr <- results(dds, contrast = c("Condition","Infected_NEG_24hr", "Infected_ILRUN_24hr"))


# Less LT07 Convert results to dataframe:
CoV2genesNEG6hr_LT07 <- results(dds_LT07, contrast = c("Condition", "Uninfected_NEG_6hr", "Infected_NEG_6hr"))
CoV2genesNEG24hr_LT07 <- results(dds_LT07, contrast = c("Condition", "Uninfected_NEG_24hr", "Infected_NEG_24hr"))
CoV2genesILRUN6hr_LT07 <- results(dds_LT07, contrast = c("Condition", "Uninfected_ILRUN_6hr", "Infected_ILRUN_6hr"))
CoV2genesILRUN24hr_LT07 <- results(dds_LT07, contrast = c("Condition", "Uninfected_ILRUN_24hr", "Infected_ILRUN_24hr"))
ILRUNgenesUninf6hr_LT07 <- results(dds_LT07, contrast = c("Condition","Uninfected_NEG_6hr", "Uninfected_ILRUN_6hr"))
ILRUNgenesUninf24hr_LT07 <- results(dds_LT07, contrast = c("Condition","Uninfected_NEG_24hr", "Uninfected_ILRUN_24hr"))
ILRUNgenesInf6hr_LT07 <- results(dds_LT07, contrast = c("Condition","Infected_NEG_6hr", "Infected_ILRUN_6hr"))
ILRUNgenesInf24hr_LT07 <- results(dds_LT07, contrast = c("Condition","Infected_NEG_24hr", "Infected_ILRUN_24hr"))

# Set adjusted p-value significance (padj) threshold:
alpha <- c( 0.05 )
# Set log2FoldChange threshold:
# As it's log2, > 1 is actually equal to 2-fold change or above.
beta <- c( 1 )
# Set baseMean threshold:
gamma <- c( 10 )

########################## CoV2genesNEG6hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesNEG6hr_LT07 <- CoV2genesNEG6hr_LT07[ which( CoV2genesNEG6hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesNEG6hr_LT07 <- sigCoV2genesNEG6hr_LT07[ which( abs(sigCoV2genesNEG6hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesNEG6hr_LT07 <- sigCoV2genesNEG6hr_LT07[ which(sigCoV2genesNEG6hr_LT07$baseMean > gamma), ]
write.csv(sigCoV2genesNEG6hr_LT07, file="results/sigCoV2genesNEG6hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesNEG6hr_LT07_hits <- rownames(sigCoV2genesNEG6hr_LT07)
length(CoV2genesNEG6hr_LT07_hits)
write.csv(norm_counts[CoV2genesNEG6hr_LT07_hits, ], file="results/DESeq2_sig_CoV2genesNEG6hr_LT07_normalized_counts.csv")
read_csv("results/sigCoV2genesNEG6hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>% 
  filter(gene != "Unknown") %>% 
  write.csv("results/CoV2genesNEG6hr_LT07_david.csv")
#making a PCA  
#read in the metadata
CoV2genesNEG6hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07.csv") %>% 
  filter(siRNA == "NEG") %>% 
  filter(timepoint == "6hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesNEG6hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesNEG6hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesNEG6hr_expression_scaled_genes <- CoV2genesNEG6hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
CoV2genesNEG6hr_pca_genes <- prcomp(CoV2genesNEG6hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_CoV2genesNEG6hr_pca_genes <- CoV2genesNEG6hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(CoV2genesNEG6hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_CoV2genesNEG6hr <- ggplot(PCA_data_CoV2genesNEG6hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Infection), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis CoV2genesNEG6hr")

ggsave(filename = "results/PCA_plot_CoV2genesNEG6hr.png", plot = plot_CoV2genesNEG6hr, width = 20, height = 20, dpi = 300, units = "cm")


########################## CoV2genesNEG24hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesNEG24hr_LT07 <- CoV2genesNEG24hr_LT07[ which( CoV2genesNEG24hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesNEG24hr_LT07 <- sigCoV2genesNEG24hr_LT07[ which( abs(sigCoV2genesNEG24hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesNEG24hr_LT07 <- sigCoV2genesNEG24hr_LT07[ which(sigCoV2genesNEG24hr_LT07$baseMean > gamma), ]
write.csv(sigCoV2genesNEG24hr_LT07, file="results/sigCoV2genesNEG24hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesNEG24hr_LT07_hits <- rownames(sigCoV2genesNEG24hr_LT07)
length(CoV2genesNEG24hr_LT07_hits)
write.csv(norm_counts[CoV2genesNEG24hr_LT07_hits, ], file="results/DESeq2_sig_CoV2genesNEG24hr_LT07_normalized_counts.csv")
read_csv("results/sigCoV2genesNEG24hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>% 
  filter(gene != "Unknown") %>%
  write.csv("results/CoV2genesNEG24hr_LT07_david.csv")
#making a PCA  
#read in the metadata
CoV2genesNEG24hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07.csv") %>% 
  filter(siRNA == "NEG") %>% 
  filter(timepoint == "24hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesNEG24hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesNEG24hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesNEG24hr_expression_scaled_genes <- CoV2genesNEG24hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
CoV2genesNEG24hr_pca_genes <- prcomp(CoV2genesNEG24hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_CoV2genesNEG24hr_pca_genes <- CoV2genesNEG24hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(CoV2genesNEG24hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_CoV2genesNEG24hr <- ggplot(PCA_data_CoV2genesNEG24hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Infection), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis CoV2genesNEG24hr")

ggsave(filename = "results/PCA_plot_CoV2genesNEG24hr.png", plot = plot_CoV2genesNEG24hr, width = 20, height = 20, dpi = 300, units = "cm")


########################## CoV2genesILRUN6hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesILRUN6hr_LT07 <- CoV2genesILRUN6hr_LT07[ which( CoV2genesILRUN6hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesILRUN6hr_LT07 <- sigCoV2genesILRUN6hr_LT07[ which( abs(sigCoV2genesILRUN6hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesILRUN6hr_LT07 <- sigCoV2genesILRUN6hr_LT07[ which(sigCoV2genesILRUN6hr_LT07$baseMean > gamma), ]
write.csv(sigCoV2genesILRUN6hr_LT07, file="results/sigCoV2genesILRUN6hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesILRUN6hr_LT07_hits <- rownames(sigCoV2genesILRUN6hr_LT07)
length(CoV2genesILRUN6hr_LT07_hits)
write.csv(norm_counts[CoV2genesILRUN6hr_LT07_hits, ], file="results/DESeq2_sig_CoV2genesILRUN6hr_LT07_normalized_counts.csv")
read_csv("results/sigCoV2genesILRUN6hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>% 
  filter(gene != "Unknown") %>% 
  write.csv("results/CoV2genesILRUN6hr_LT07_david.csv")
#making a PCA  
#read in the metadata
CoV2genesILRUN6hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07.csv") %>% 
  filter(siRNA == "ILRUN") %>% 
  filter(timepoint == "6hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesILRUN6hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesILRUN6hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesILRUN6hr_expression_scaled_genes <- CoV2genesILRUN6hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
CoV2genesILRUN6hr_pca_genes <- prcomp(CoV2genesILRUN6hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_CoV2genesILRUN6hr_pca_genes <- CoV2genesILRUN6hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(CoV2genesILRUN6hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_CoV2genesILRUN6hr <- ggplot(PCA_data_CoV2genesILRUN6hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Infection), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis CoV2genesILRUN6hr")

ggsave(filename = "results/PCA_plot_CoV2genesILRUN6hr.png", plot = plot_CoV2genesILRUN6hr, width = 20, height = 20, dpi = 300, units = "cm")


########################## CoV2genesILRUN24hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesILRUN24hr_LT07 <- CoV2genesILRUN24hr_LT07[ which( CoV2genesILRUN24hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesILRUN24hr_LT07 <- sigCoV2genesILRUN24hr_LT07[ which( abs(sigCoV2genesILRUN24hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesILRUN24hr_LT07 <- sigCoV2genesILRUN24hr_LT07[ which(sigCoV2genesILRUN24hr_LT07$baseMean > gamma), ]
write.csv(sigCoV2genesILRUN24hr_LT07, file="results/sigCoV2genesILRUN24hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesILRUN24hr_LT07_hits <- rownames(sigCoV2genesILRUN24hr_LT07)
length(CoV2genesILRUN24hr_LT07_hits)
write.csv(norm_counts[CoV2genesILRUN24hr_LT07_hits, ], file="results/DESeq2_sig_CoV2genesILRUN24hr_LT07_normalized_counts.csv")
read_csv("results/sigCoV2genesILRUN24hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>% 
  filter(gene != "Unknown") %>% 
  write.csv("results/CoV2genesILRUN24hr_LT07_david.csv")
#making a PCA  
#read in the metadata
CoV2genesILRUN24hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07.csv") %>% 
  filter(siRNA == "ILRUN") %>% 
  filter(timepoint == "24hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesILRUN24hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesILRUN24hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesILRUN24hr_expression_scaled_genes <- CoV2genesILRUN24hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
CoV2genesILRUN24hr_pca_genes <- prcomp(CoV2genesILRUN24hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_CoV2genesILRUN24hr_pca_genes <- CoV2genesILRUN24hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(CoV2genesILRUN24hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_CoV2genesILRUN24hr <- ggplot(PCA_data_CoV2genesILRUN24hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Infection), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis CoV2genesILRUN24hr")

ggsave(filename = "results/PCA_plot_CoV2genesILRUN24hr.png", plot = plot_CoV2genesILRUN24hr, width = 20, height = 20, dpi = 300, units = "cm")

########################## ILRUNgenesUninf6hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigILRUNgenesUninf6hr_LT07<- ILRUNgenesUninf6hr_LT07[ which( ILRUNgenesUninf6hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesUninf6hr_LT07 <- sigILRUNgenesUninf6hr_LT07[ which( abs(sigILRUNgenesUninf6hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesUninf6hr_LT07 <- sigILRUNgenesUninf6hr_LT07[ which(sigILRUNgenesUninf6hr_LT07$baseMean > gamma), ]
write.csv(sigILRUNgenesUninf6hr_LT07, file="results/sigILRUNgenesUninf6hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesUninf6hr_LT07_hits <- rownames(sigILRUNgenesUninf6hr_LT07)
length(ILRUNgenesUninf6hr_LT07_hits)
write.csv(norm_counts[ILRUNgenesUninf6hr_LT07_hits, ], file="results/DESeq2_sig_ILRUNgenesUninf6hr_LT07_normalized_counts.csv")
read_csv("results/sigILRUNgenesUninf6hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>%
  filter(gene != "Unknown") %>%
  write.csv("results/ILRUNgenesUninf6hr_LT07_david.csv")

########################## ILRUNgenesUninf24hr #####################################################################

# Filter ILRUNgenesUninf24hr:
# 'Which' provides positions of 'TRUE' values.
sigILRUNgenesUninf24hr_LT07<- ILRUNgenesUninf24hr_LT07[ which( ILRUNgenesUninf24hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesUninf24hr_LT07 <- sigILRUNgenesUninf24hr_LT07[ which( abs(sigILRUNgenesUninf24hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesUninf24hr_LT07 <- sigILRUNgenesUninf24hr_LT07[ which(sigILRUNgenesUninf24hr_LT07$baseMean > gamma), ]
write.csv(sigILRUNgenesUninf24hr_LT07, file="results/sigILRUNgenesUninf24hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesUninf24hr_LT07_hits <- rownames(sigILRUNgenesUninf24hr_LT07)
length(ILRUNgenesUninf24hr_LT07_hits)
write.csv(norm_counts[ILRUNgenesUninf24hr_LT07_hits, ], file="results/DESeq2_sig_ILRUNgenesUninf24hr_LT07_normalized_counts.csv")
read_csv("results/sigILRUNgenesUninf24hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>%
  filter(gene != "Unknown") %>%
  write.csv("results/ILRUNgenesUninf24hr_LT07_david.csv")

########################## ILRUNgenesInf6hr ###################################################################### 
##'Which' provides positions of 'TRUE' values.
sigILRUNgenesInf6hr_LT07<- ILRUNgenesInf6hr_LT07[ which( ILRUNgenesInf6hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesInf6hr_LT07 <- sigILRUNgenesInf6hr_LT07[ which( abs(sigILRUNgenesInf6hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesInf6hr_LT07 <- sigILRUNgenesInf6hr_LT07[ which(sigILRUNgenesInf6hr_LT07$baseMean > gamma), ]
write.csv(sigILRUNgenesInf6hr_LT07, file="results/sigILRUNgenesInf6hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesInf6hr_LT07_hits <- rownames(sigILRUNgenesInf6hr_LT07)
length(ILRUNgenesInf6hr_LT07_hits)
write.csv(norm_counts[ILRUNgenesInf6hr_LT07_hits, ], file="results/DESeq2_sig_ILRUNgenesInf6hr_LT07_normalized_counts.csv")
read_csv("results/sigILRUNgenesInf6hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>%
  filter(gene != "Unknown") %>%
  write.csv("results/ILRUNgenesInf6hr_LT07_david.csv")

########################## ILRUNgenesInf24hr ###################################################################### 
# 'Which' provides positions of 'TRUE' values.
sigILRUNgenesInf24hr_LT07<- ILRUNgenesInf24hr_LT07[ which( ILRUNgenesInf24hr_LT07$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesInf24hr_LT07 <- sigILRUNgenesInf24hr_LT07[ which( abs(sigILRUNgenesInf24hr_LT07$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesInf24hr_LT07 <- sigILRUNgenesInf24hr_LT07[ which(sigILRUNgenesInf24hr_LT07$baseMean > gamma), ]
write.csv(sigILRUNgenesInf24hr_LT07, file="results/sigILRUNgenesInf24hr_LT07.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesInf24hr_LT07_hits <- rownames(sigILRUNgenesInf24hr_LT07)
length(ILRUNgenesInf24hr_LT07_hits)
write.csv(norm_counts[ILRUNgenesInf24hr_LT07_hits, ], file="results/DESeq2_sig_ILRUNgenesInf24hr_LT07_normalized_counts.csv")
read_csv("results/sigILRUNgenesInf24hr_LT07.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(col =gene, into = c("name", "XLOC", "gene"), sep ="_") %>% 
  select(gene) %>%
  filter(gene != "Unknown") %>%
  write.csv("results/ILRUNgenesInf24hr_LT07_david.csv")






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
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
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
  geom_point(aes(shape = siRNA, color = Infection), size = 6)+
  geom_text(aes(label = timepoint), size = 4) +
  labs(title = "Principle Component Analysis less LT07")
  
ggsave(filename = "results/PCA_lessLT07.png", plot = text_plot_lessLT07, width = 12, height = 10, dpi = 300, units = "cm")




#####Repeat analysis for dataset excluding LT07 and excluding lane duplicates for plot readability###############
#read in the metadata
lessLT07_lessLane_metadata <- read_csv("data/ilrun_metadata_siRNA.csv") %>% 
  mutate(Name = str_extract(Sample, "^....")) %>%
  filter(Sample != "LT07_L001" |Sample != "LT07_L002") %>%
  filter(Lane != "L002")
  
#join with the counts data 
lessLT07_lessLane_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts.csv") %>%
  rename(gene_locus = X1) %>%
  select(-LT07_L001, -LT07_L002) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  left_join(lessLT07_lessLane_metadata, by = "Sample")

#scale gene expression for all samples 
lessLT07_lessLane_scaled_genes <- lessLT07_lessLane_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
lessLT07_Lesslane_pca_genes <- prcomp(lessLT07_lessLane_scaled_genes)

#tidy data frame for plotting
PCA_data_lessLT07_Lesslane <- lessLT07_Lesslane_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(all_metadata, by = "Sample") %>%
  spread(PC, expression)

text_plot_lessLT07_Lesslane <- ggplot(PCA_data_lessLT07_Lesslane, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = siRNA, color = Infection), size = 6)+
  geom_text(aes(label = timepoint), size = 4) +
  labs(title = "Principle Component Analysis less LT07 less Lanes all normalized counts")

ggsave(filename = "results/PCA_lessLT07_LessLane.png", plot = text_plot_lessLT07_Lesslane, width = 12, height = 10, dpi = 300, units = "cm")




######### Expression of ISGs in ILRUN data ######

write.csv(lessLT07_expression, file="results/DESeq2_ilrun_siRNA_lessLT07_expression.csv")


genes_of_interest <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_000014_ISG15" | gene_locus == "XLOC_001047_IL6R" | gene_locus == "XLOC_001145_FCER1G" |
           gene_locus == "XLOC_021735_IFNAR1" | gene_locus == "XLOC_036412_ACE2")
  

plot_genes_of_interest <- ggplot(genes_of_interest, aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot() +
  facet_wrap(~gene_locus)

ISG15 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_000014_ISG15") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "ISG15", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

IL6R <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_001047_IL6R") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "IL6R", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

FCER1G <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_001145_FCER1G") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "FCER1G", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

ACE2 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_036412_ACE2") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "ACE2", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

IFNAR1 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_021735_IFNAR1") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "IFNAR1", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

VIM <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_003358_VIM") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "VIM", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

GBP1 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_002271_GBP1") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "GBP1", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))


FOLR2 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_005291_FOLR2") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "FOLR2", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

IL15RA <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_004013_IL15RA") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "IL15RA", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

library(cowplot)

# extract the legend from one of the plots
legend <- get_legend(ISG15) 

plot <- plot_grid(ISG15 + theme(legend.position="none"),
                  IL6R + theme(legend.position="none"),
                  FCER1G + theme(legend.position="none"),
                  ACE2 + theme(legend.position="none"),
                  IFNAR1 + theme(legend.position="none"),
                  VIM + theme(legend.position="none"), 
                  GBP1 +theme(legend.position="none"),
                  FOLR2 +theme(legend.position="none"),
                  IL15RA +theme(legend.position="none"))

genes_of_interest <- plot_grid(plot, legend, rel_widths = c(3, .3))

ggsave(filename = "results/genes_of_interest.png", plot = genes_of_interest, width = 25, height = 20, dpi = 300, units = "cm")



FBL <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_017255_FBL") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "FBL", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))



RPL13 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_012366_RPL13") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "RPL13", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))


RPL7A <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_034787_RPL7A") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "RPL7A", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))


RPL4 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_011308_RPL4") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "RPL4", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

RPL13A <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_016495_RPL13A") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "RPL13A", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))


SNORD32A <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_016496_SNORD32A") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "SNORD32A", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

SNORD34 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_016498_SNORD34") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "SNORD34", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

SNORD35A <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_016499_SNORD35A") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "SNORD35A", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

SNORD36A <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_034790_SNORD36A") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "SNORD36A", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

SNORD36C <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_034791_SNORD36C") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "SNORD36C", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

SNORD38B <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_000458_SNORD38B") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "SNORD38B", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

PDZD4 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_037006_PDZD4") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "PDZD4", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

DPP4 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_020331_DPP4") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "DPP4", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))


TRIM31 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_028190_TRIM31") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "TRIM31", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

ggsave(filename = "results/TRIM31.png", plot = TRIM31, width = 25, height = 20, dpi = 300, units = "cm")



LGI4 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_017175_LGI4") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "LGI4", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

ggsave(filename = "results/LGI4.png", plot = LGI4, width = 25, height = 20, dpi = 300, units = "cm")

PDZD4 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_037006_PDZD4") %>% 
  ggplot(aes(x= siRNA, y = expression)) +
  geom_boxplot()+
  labs(title = "PDZD4", 
       x = "")
