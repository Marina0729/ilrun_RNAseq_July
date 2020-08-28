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
library(EnhancedVolcano)
library(viridis)

rm(list=ls())

#splice the count data 
read.csv("data/all_samples_rev-intersect.csv") %>% 
  select(-LT41_L001_rev.intersect) %>% 
  select(-LT41_L002_rev.intersect) %>% 
  select(-LT42_L001_rev.intersect) %>% 
  select(-LT42_L002_rev.intersect) %>% 
  select(-LT07_L001_rev.intersect) %>% 
  select(-LT07_L002_rev.intersect) %>%
  as_tibble() %>% 
  write.csv("data/ilrun_counts_siRNA_LT07_rev.csv",  row.names=FALSE)

df <- read.csv("data/ilrun_counts_siRNA_LT07_rev.csv", header=TRUE, row.names=1)

dim(df)

#splice the metadata
ilrun_metadata_siRNA_rev_LT07 <- read.csv("data/ilrun_metadata_siRNA_LT07.csv") %>%
  mutate(HTSeq_rev = "rev.intersect" ) %>% 
  unite("Sample", Sample, HTSeq_rev, remove = FALSE) %>% 
  write.csv("data/ilrun_metadata_siRNA_LT07_rev.csv")

md <- read.csv("data/ilrun_metadata_siRNA_LT07_rev.csv", header=TRUE, row.names=2)


dim(md)

#merge count and metadata file
all(rownames(md) == colnames (df))

# Calculate counts per million.
# Filter rows: at least 3 samples must have at least 1 cpm.
# Retain rows according to 'keep' and nothing on columns.
cpm <- apply(df, 2, function(x) (x/sum(x))*1000000)
keep <- rowSums( cpm >= 1 ) >=4
df_filtered <- df[ keep, ]

dim(df_filtered)


# Construct a SummarizedExperiment object:
dds <- DESeqDataSetFromMatrix(
  countData = df_filtered,
  colData = md,
  design = ~ Condition + Lane) # ~ is representative of 'by', i.e. compare by condition


# Perform DE testing:
dds <- DESeq(dds)

# Output normalized counts:
norm_counts <- counts (dds, normalized=TRUE)
write.csv(norm_counts, file="results/DESeq2_ilrun_siRNA_LT07_rev_all_normalized_counts.csv")


# Convert results to dataframe:
CoV2genesNEG6hr <- results(dds, contrast = c("Condition", "Infected_NEG_6hr", "Uninfected_NEG_6hr"))
CoV2genesNEG24hr <- results(dds, contrast = c("Condition", "Infected_NEG_24hr", "Uninfected_NEG_24hr"))
CoV2genesILRUN6hr <- results(dds, contrast = c("Condition", "Infected_ILRUN_6hr", "Uninfected_ILRUN_6hr"))
CoV2genesILRUN24hr <- results(dds, contrast = c("Condition", "Infected_ILRUN_24hr", "Uninfected_ILRUN_24hr"))

ILRUNgenesUninf6hr <- results(dds, contrast = c("Condition", "Uninfected_ILRUN_6hr","Uninfected_NEG_6hr"))
ILRUNgenesUninf24hr <- results(dds, contrast = c("Condition", "Uninfected_ILRUN_24hr","Uninfected_NEG_24hr"))
ILRUNgenesInf6hr <- results(dds, contrast = c("Condition", "Infected_ILRUN_6hr", "Infected_NEG_6hr"))
ILRUNgenesInf24hr <- results(dds, contrast = c("Condition", "Infected_ILRUN_24hr", "Infected_NEG_24hr"))


# Set adjusted p-value significance (padj) threshold:
alpha <- c( 0.05 )
# Set log2FoldChange threshold:
# As it's log2, > 1 is actually equal to 2-fold change or above.
beta <- c( 0.75 )
# Set baseMean threshold:
gamma <- c( 5 )

########################## CoV2genesNEG6hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesNEG6hr <- CoV2genesNEG6hr[ which( CoV2genesNEG6hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesNEG6hr <- sigCoV2genesNEG6hr[ which( abs(sigCoV2genesNEG6hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesNEG6hr <- sigCoV2genesNEG6hr[ which(sigCoV2genesNEG6hr$baseMean > gamma), ]
write.csv(sigCoV2genesNEG6hr, file="results/sigCoV2genesNEG6hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesNEG6hr_hits <- rownames(sigCoV2genesNEG6hr)
length(CoV2genesNEG6hr_hits)
write.csv(norm_counts[CoV2genesNEG6hr_hits, ], file="results//DESeq2_sig_CoV2genesNEG6hr_LT07_rev_normalized_counts")
read_csv("results/sigCoV2genesNEG6hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/CoV2genesNEG6hr_LT07_rev_david.csv")

#making a PCA  
#read in the metadata
CoV2genesNEG6hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(siRNA == "NEG") %>% 
  filter(timepoint == "6hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesNEG6hr_expression <- read_csv("results//DESeq2_ilrun_siRNA_LT07_rev_all_normalized_counts.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesNEG6hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesNEG6hr_expression_scaled_genes <- CoV2genesNEG6hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
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

forplotting_CoV2genesNEG6hr <- CoV2genesNEG6hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>%
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")
  
volcano_plot_CoV2genesNEG6hr <- EnhancedVolcano(forplotting_CoV2genesNEG6hr,
                lab = rownames(forplotting_CoV2genesNEG6hr),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-3, 3),
                title = 'CoV2genes siNEG 6hr',
                pCutoff = 0.05,
                FCcutoff = 0.75,
                pointSize = 3.0,
                labSize = 3.0)

ggsave(filename = "results/volcano_plot_CoV2genesNEG6hr.png", plot = volcano_plot_CoV2genesNEG6hr, width = 20, height = 20, dpi = 300, units = "cm")


########################## CoV2genesNEG24hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesNEG24hr <- CoV2genesNEG24hr[ which( CoV2genesNEG24hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesNEG24hr <- sigCoV2genesNEG24hr[ which( abs(sigCoV2genesNEG24hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesNEG24hr <- sigCoV2genesNEG24hr[ which(sigCoV2genesNEG24hr$baseMean > gamma), ]
write.csv(sigCoV2genesNEG24hr, file="results/sigCoV2genesNEG24hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesNEG24hr_hits <- rownames(sigCoV2genesNEG24hr)
length(CoV2genesNEG24hr_hits)
write.csv(norm_counts[CoV2genesNEG24hr_hits, ], file="results/DESeq2_sig_CoV2genesNEG24hr_LT07_rev_normalized_counts")
read_csv("results/sigCoV2genesNEG24hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/CoV2genesNEG24hr_LT07_david.csv")
#making a PCA  
#read in the metadata
CoV2genesNEG24hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(siRNA == "NEG") %>% 
  filter(timepoint == "24hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesNEG24hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_LT07_rev_all_normalized_counts.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesNEG24hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesNEG24hr_expression_scaled_genes <- CoV2genesNEG24hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
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

forplotting_CoV2genesNEG24hr <- CoV2genesNEG24hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>%
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

volcano_plot_CoV2genesNEG24hr <- EnhancedVolcano(forplotting_CoV2genesNEG24hr,
                                                lab = rownames(forplotting_CoV2genesNEG24hr),
                                                x = 'log2FoldChange',
                                                y = 'pvalue',
                                                xlim = c(-3, 3),
                                                title = 'CoV2genes siNEG 24hr',
                                                pCutoff = 0.05,
                                                FCcutoff = 0.75,
                                                pointSize = 3.0,
                                                labSize = 3.0)

ggsave(filename = "results/volcano_plot_CoV2genesNEG624hr.png", plot = volcano_plot_CoV2genesNEG24hr, width = 20, height = 20, dpi = 300, units = "cm")


########################## CoV2genesILRUN6hr #####################################################################
# 'Which' provides positions of 'TRUE' values
sigCoV2genesILRUN6hr <- CoV2genesILRUN6hr[ which( CoV2genesILRUN6hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesILRUN6hr <- sigCoV2genesILRUN6hr[ which( abs(sigCoV2genesILRUN6hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesILRUN6hr <- sigCoV2genesILRUN6hr[ which(sigCoV2genesILRUN6hr$baseMean > gamma), ]
write.csv(sigCoV2genesILRUN6hr, file="results/sigCoV2genesILRUN6hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesILRUN6hr_hits <- rownames(sigCoV2genesILRUN6hr)
length(CoV2genesILRUN6hr_hits)
write.csv(norm_counts[CoV2genesILRUN6hr_hits, ], file="results/DESeq2_sig_CoV2genesILRUN6hr_LT07_rev_normalized_counts.csv")
read_csv("results/sigCoV2genesILRUN6hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/CoV2genesILRUN6hr_LT07_rev_david.csv")
#making a PCA  
#read in the metadata
CoV2genesILRUN6hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(siRNA == "ILRUN") %>% 
  filter(timepoint == "6hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesILRUN6hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07_rev.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesILRUN6hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesILRUN6hr_expression_scaled_genes <- CoV2genesILRUN6hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
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

forplotting_CoV2genesILRUN6hr <- CoV2genesILRUN6hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>%
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

volcano_plot_CoV2genesILRUN6hr <- EnhancedVolcano(forplotting_CoV2genesNEG6hr,
                                                 lab = rownames(forplotting_CoV2genesILRUN6hr),
                                                 x = 'log2FoldChange',
                                                 y = 'pvalue',
                                                 xlim = c(-3, 3),
                                                 title = 'CoV2genes siILRUN 6hr',
                                                 pCutoff = 0.05,
                                                 FCcutoff = 0.75,
                                                 pointSize = 3.0,
                                                 labSize = 3.0)

ggsave(filename = "results/volcano_plot_CoV2genesILRUN6hr.png", plot = volcano_plot_CoV2genesILRUN6hr, width = 20, height = 20, dpi = 300, units = "cm")

########################## CoV2genesILRUN24hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigCoV2genesILRUN24hr <- CoV2genesILRUN24hr[ which( CoV2genesILRUN24hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigCoV2genesILRUN24hr <- sigCoV2genesILRUN24hr[ which( abs(sigCoV2genesILRUN24hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigCoV2genesILRUN24hr <- sigCoV2genesILRUN24hr[ which(sigCoV2genesILRUN24hr$baseMean > gamma), ]
write.csv(sigCoV2genesILRUN24hr, file="results/sigCoV2genesILRUN24hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
CoV2genesILRUN24hr_hits <- rownames(sigCoV2genesILRUN24hr)
length(CoV2genesILRUN24hr_hits)
write.csv(norm_counts[CoV2genesILRUN24hr_hits, ], file="results/DESeq2_sig_CoV2genesILRUN24hr_LT07_rev_normalized_counts.csv")
read_csv("results/sigCoV2genesILRUN24hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/CoV2genesILRUN24hr_LT07_david_rev.csv")
#making a PCA  
#read in the metadata
CoV2genesILRUN24hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(siRNA == "ILRUN") %>% 
  filter(timepoint == "24hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
CoV2genesILRUN24hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07_rev.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(CoV2genesILRUN24hr_metadata, by = "Sample") 

#scale gene expression for all samples 
CoV2genesILRUN24hr_expression_scaled_genes <- CoV2genesILRUN24hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
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

forplotting_CoV2genesILRUN24hr <- CoV2genesILRUN24hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

volcano_plot_CoV2genesILRUN24hr <- EnhancedVolcano(forplotting_CoV2genesILRUN24hr,
                                                 lab = rownames(forplotting_CoV2genesILRUN24hr),
                                                 x = 'log2FoldChange',
                                                 y = 'pvalue',
                                                 xlim = c(-3, 3),
                                                 title = 'CoV2genes siILRUN 24hr',
                                                 pCutoff = 0.05,
                                                 FCcutoff = 0.75,
                                                 pointSize = 3.0,
                                                 labSize = 3.0)

ggsave(filename = "results/volcano_plot_CoV2genesILRUN24hr.png", plot = volcano_plot_CoV2genesILRUN24hr, width = 20, height = 20, dpi = 300, units = "cm")

########################## ILRUNgenesUninf6hr #####################################################################
# 'Which' provides positions of 'TRUE' values.
sigILRUNgenesUninf6hr<- ILRUNgenesUninf6hr[ which( ILRUNgenesUninf6hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesUninf6hr <- sigILRUNgenesUninf6hr[ which( abs(sigILRUNgenesUninf6hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesUninf6hr <- sigILRUNgenesUninf6hr[ which(sigILRUNgenesUninf6hr$baseMean > gamma), ]
write.csv(sigILRUNgenesUninf6hr, file="results/sigILRUNgenesUninf6hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesUninf6hr_hits <- rownames(sigILRUNgenesUninf6hr)
length(ILRUNgenesUninf6hr_hits)
write.csv(norm_counts[ILRUNgenesUninf6hr_hits, ], file="results/DESeq2_sig_ILRUNgenesUninf6hr_LT07_rev_normalized_counts.csv")
read_csv("results/sigILRUNgenesUninf6hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/ILRUNgenesUninf6hr_LT07_rev_david.csv")

######################## making a PCA ##########################################
#read in the metadata
ILRUNgenesUninf6hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(Infection == "Uninfected") %>% 
  filter(timepoint == "6hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
ILRUNgenesUninf6hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07_rev.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(ILRUNgenesUninf6hr_metadata, by = "Sample") 

#scale gene expression for all samples 
ILRUNgenesUninf6hr_expression_scaled_genes <- ILRUNgenesUninf6hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
ILRUNgenesUninf6hr_pca_genes <- prcomp(ILRUNgenesUninf6hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_ILRUNgenesUninf6hr_pca_genes <- ILRUNgenesUninf6hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(ILRUNgenesUninf6hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_ILRUNgenesUninf6hr <- ggplot(PCA_data_ILRUNgenesUninf6hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = siRNA), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis ILRUNgenesUninf6hr")

ggsave(filename = "results/PCA_plot_ILRUNgenesUninf6hr.png", plot = plot_ILRUNgenesUninf6hr, width = 20, height = 20, dpi = 300, units = "cm")

######################## Volcano plots ###################################

forplotting_ILRUNgenesUninf6hr <- ILRUNgenesUninf6hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  mutate(gene = str_replace(gene, "ACKR3", "CXCR7")) %>%
  mutate(gene = str_replace(gene, "CHUK","IkBKA")) %>%
  mutate(gene = str_replace(gene, "TNFSF10", "TRAIL")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

write.csv(forplotting_ILRUNgenesUninf6hr, "results/ILRUNgenesUninf6hr_forplotting.csv")

##show only significant immune genes ###
tibble_ILRUNgenesUninf6hr <- ILRUNgenesUninf6hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  mutate(gene = str_replace(gene, "ACKR3", "CXCR7")) %>%
  mutate(gene = str_replace(gene, "CHUK","IkBKA")) %>%
  mutate(gene = str_replace(gene, "TNFSF10", "TRAIL"))


immune_genes_Uninf6hr <- c('ISG15','IL6R','TLR7', 'STAT4', 'SOCS2','IL1B','IL15RA', 'AKT1', 'CD14', 'TLR2', 'IFIH1', 'OAS1',
                  'JAK1', 'JAK2', 'CDK2', 'NFKB1', 'IL12A', 'IL8', 'TAP1', 'HLA-B', 'HLA-DMA', 'TRAIL', 'IL2RG', 'IL1R',
                  'TGFB1', 'CCR1', 'CCR6', 'IL8', 'IFIT3', 'IL18', 'IFIT3', 'CXCR7', 'GDF11', 'LIFR', 'MMP14', 'THBS1', 'VAV3', 'PRKCQ',
                  'IkBKA', 'ILRUN', 'VIM', 'FCGR2') %>%
  as_tibble() %>% 
  rename(gene = value) %>% 
  left_join(tibble_ILRUNgenesUninf6hr , by = "gene") %>% 
  filter(padj <0.05 ) %>% 
  filter(log2FoldChange >0.75 | log2FoldChange < -0.75) 
  

#########################################

volcano_plot_ILRUNgenesUninf6hr <- EnhancedVolcano(forplotting_ILRUNgenesUninf6hr,
                                                 lab = rownames(forplotting_ILRUNgenesUninf6hr),
                                                 x = 'log2FoldChange',
                                                 y = 'padj',
                                                 selectLab = c(immune_genes_Uninf6hr$gene),
                                                 xlim = c(-3, 3),
                                                 title = 'ILRUN genes 6hr',
                                                 pCutoff = 0.05,
                                                 FCcutoff = 0.5,
                                                 pointSize = 4.0,
                                                 labSize = 5.0,
                                                 labCol = 'black',
                                                 labFace = 'bold',
                                                 boxedLabels = TRUE,
                                                 colAlpha = 4/5,
                                                 legendPosition = 'right',
                                                 legendLabSize = 14,
                                                 legendIconSize = 4.0,
                                                 drawConnectors = TRUE,
                                                 widthConnectors = 1.0,
                                                 colConnectors = 'black')


ggsave(filename = "results/volcano_plot_ILRUNgenesUninf6hr.png", plot = volcano_plot_ILRUNgenesUninf6hr, width = 30, height = 25, dpi = 300, units = "cm")

########################## ILRUNgenesUninf24hr #####################################################################

# Filter ILRUNgenesUninf24hr:
# 'Which' provides positions of 'TRUE' values.
sigILRUNgenesUninf24hr<- ILRUNgenesUninf24hr[ which( ILRUNgenesUninf24hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesUninf24hr <- sigILRUNgenesUninf24hr[ which( abs(sigILRUNgenesUninf24hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesUninf24hr <- sigILRUNgenesUninf24hr[ which(sigILRUNgenesUninf24hr$baseMean > gamma), ]
write.csv(sigILRUNgenesUninf24hr, file="results/sigILRUNgenesUninf24hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesUninf24hr_hits <- rownames(sigILRUNgenesUninf24hr)
length(ILRUNgenesUninf24hr_hits)
write.csv(norm_counts[ILRUNgenesUninf24hr_hits, ], file="results/DESeq2_sig_ILRUNgenesUninf24hr_LT07_rev_normalized_counts.csv")

david <- read_csv("results/sigILRUNgenesUninf24hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  filter(str_detect(gene, "LOC") == FALSE) %>% 
  write.csv("results/ILRUNgenesUninf24hr_LT07_rev_david.csv")

#making a PCA  
#read in the metadata
ILRUNgenesUninf24hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(Infection == "Uninfected") %>% 
  filter(timepoint == "24hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
ILRUNgenesUninf24hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07_rev.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(ILRUNgenesUninf24hr_metadata, by = "Sample") 

#scale gene expression for all samples 
ILRUNgenesUninf24hr_expression_scaled_genes <- ILRUNgenesUninf24hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
ILRUNgenesUninf24hr_pca_genes <- prcomp(ILRUNgenesUninf24hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_ILRUNgenesUninf24hr_pca_genes <- ILRUNgenesUninf24hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(ILRUNgenesUninf24hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_ILRUNgenesUninf24hr <- ggplot(PCA_data_ILRUNgenesUninf24hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = siRNA), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis ILRUNgenesUninf24hr")

ggsave(filename = "results/PCA_plot_ILRUNgenesUninf24hr.png", plot = plot_ILRUNgenesUninf24hr, width = 20, height = 20, dpi = 300, units = "cm")

######################## Volcano plots ###################################

forplotting_ILRUNgenesUninf24hr <- ILRUNgenesUninf24hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  mutate(gene = str_replace(gene, "ACKR3", "CXCR7")) %>%
  mutate(gene = str_replace(gene, "CHUK","IkBKA")) %>%
  mutate(gene = str_replace(gene, "TNFSF10", "TRAIL")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

write.csv(forplotting_ILRUNgenesUninf24hr, "results/ILRUNgenesUninf24hr_forplotting.csv")

##show only significant immune genes ###
tibble_ILRUNgenesUninf24hr <- ILRUNgenesUninf24hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  mutate(gene = str_replace(gene, "ACKR3", "CXCR7")) %>%
  mutate(gene = str_replace(gene, "CHUK","IkBKA")) %>%
  mutate(gene = str_replace(gene, "TNFSF10", "TRAIL"))


immune_genes_Uninf24hr <- c('ISG15','IL6R','TLR7', 'STAT4', 'SOCS2','IL1B','IL15RA', 'AKT1', 'CD14', 'TLR2', 'IFIH1', 'OAS1',
                           'JAK1', 'JAK2', 'CDK2', 'NFKB1', 'IL12A', 'IL8', 'TAP1', 'HLA-B', 'HLA-DMA', 'TRAIL', 'IL2RG', 'IL1R',
                           'TGFB1', 'CCR1', 'CCR6', 'IL8', 'IFIT3', 'IL18', 'IFIT3', 'CXCR7', 'GDF11', 'LIFR', 'MMP14', 'THBS1', 'VAV3', 'PRKCQ',
                           'IkBKA', 'ILRUN', 'VIM', 'FCGR2') %>%
  as_tibble() %>% 
  rename(gene = value) %>% 
  left_join(tibble_ILRUNgenesUninf24hr , by = "gene") %>% 
  filter(padj <0.05 ) %>% 
  filter(log2FoldChange >0.75 | log2FoldChange < -0.75) 


#########################################

volcano_plot_ILRUNgenesUninf24hr <- EnhancedVolcano(forplotting_ILRUNgenesUninf24hr,
                                                   lab = rownames(forplotting_ILRUNgenesUninf24hr),
                                                   x = 'log2FoldChange',
                                                   y = 'padj',
                                                   selectLab = c(immune_genes_Uninf24hr$gene),
                                                   xlim = c(-3, 3),
                                                   title = 'ILRUN genes 24hr',
                                                   pCutoff = 0.05,
                                                   FCcutoff = 0.5,
                                                   pointSize = 4.0,
                                                   labSize = 5.0,
                                                   labCol = 'black',
                                                   labFace = 'bold',
                                                   boxedLabels = TRUE,
                                                   colAlpha = 4/5,
                                                   legendPosition = 'right',
                                                   legendLabSize = 14,
                                                   legendIconSize = 4.0,
                                                   drawConnectors = TRUE,
                                                   widthConnectors = 1.0,
                                                   colConnectors = 'black')

ggsave(filename = "results/volcano_plot_ILRUNgenesUninf24hr.png", plot = volcano_plot_ILRUNgenesUninf24hr, width = 30, height = 25, dpi = 300, units = "cm")


############### David Functional Annotation Clustering #############

annotation_ILRUNgenesUninf24hr <- read.csv("data/ILRUNgenes_Uninf_24hrs_david_default_output.csv") %>% 
  as_tibble() %>% 
  filter(Category == "BIOCARTA" | Category == "KEGG_PATHWAY") %>% 
  arrange(PValue) %>%
  separate(Term, c("code", "Pathway"), sep = ":") %>% 
  mutate(log10P = -log(PValue))

write.csv(annotation_ILRUNgenesUninf24hr, "results/pathway_annotation_ILRUNgenesUninf24hr.csv")

plot_ILRUNgenesUninf24hr_pathways <- ggplot(annotation_ILRUNgenesUninf24hr) + 
  geom_bar(aes(x = reorder(Pathway, log10P), y = Fold_enrichment, fill = log10P), stat = 'identity') +
  coord_flip() +
  scale_fill_viridis_c() +
  labs(title = "Pathways enriched for ILRUN regulated genes",
       y = "Fold Enrichment",
       x = "Pathway",
       fill = "-log10(p-value)")

ggsave(filename = "results/pathway_plot_ILRUNgenesUninf24hr.png", plot = plot_ILRUNgenesUninf24hr_pathways, width = 30, height = 20, dpi = 300, units = "cm")

########################## ILRUNgenesInf6hr ###################################################################### 
##'Which' provides positions of 'TRUE' values.
sigILRUNgenesInf6hr<- ILRUNgenesInf6hr[ which( ILRUNgenesInf6hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesInf6hr <- sigILRUNgenesInf6hr[ which( abs(sigILRUNgenesInf6hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesInf6hr <- sigILRUNgenesInf6hr[ which(sigILRUNgenesInf6hr$baseMean > gamma), ]
write.csv(sigILRUNgenesInf6hr, file="results/sigILRUNgenesInf6hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesInf6hr_hits <- rownames(sigILRUNgenesInf6hr)
length(ILRUNgenesInf6hr_hits)
write.csv(norm_counts[ILRUNgenesInf6hr_hits, ], file="results/DESeq2_sig_ILRUNgenesInf6hr_LT07_rev_normalized_counts.csv")
read_csv("results/sigILRUNgenesInf6hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/ILRUNgenesInf6hr_LT07_rev_david.csv")

#making a PCA  
#read in the metadata
ILRUNgenesInf6hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(Infection == "Infected") %>% 
  filter(timepoint == "6hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
ILRUNgenesInf6hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07_rev.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(ILRUNgenesInf6hr_metadata, by = "Sample") 

#scale gene expression for all samples 
ILRUNgenesInf6hr_expression_scaled_genes <- ILRUNgenesInf6hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>%
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
ILRUNgenesInf6hr_pca_genes <- prcomp(ILRUNgenesInf6hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_ILRUNgenesInf6hr_pca_genes <- ILRUNgenesInf6hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(ILRUNgenesInf6hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_ILRUNgenesInf6hr <- ggplot(PCA_data_ILRUNgenesInf6hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = siRNA), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis ILRUNgenesInf6hr")

ggsave(filename = "results/PCA_plot_ILRUNgenesInf6hr.png", plot = plot_ILRUNgenesInf6hr, width = 20, height = 20, dpi = 300, units = "cm")

forplotting_ILRUNgenesInf6hr <- ILRUNgenesInf6hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>%
  mutate(gene = str_replace(gene, "C6orf106", "ILRUN")) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")
 

volcano_plot_ILRUNgenesInf6hr <- EnhancedVolcano(forplotting_ILRUNgenesInf6hr,
                                                    lab = rownames(forplotting_ILRUNgenesInf6hr),
                                                    x = 'log2FoldChange',
                                                    y = 'pvalue',
                                                    xlim = c(-3, 3),
                                                    title = 'ILRUN genes Infected 6hr',
                                                    pCutoff = 0.05,
                                                    FCcutoff = 0.75,
                                                    pointSize = 3.0,
                                                    labSize = 3.0)

ggsave(filename = "results/volcano_plot_ILRUNgenesInf6hr.png", plot = volcano_plot_ILRUNgenesInf6hr, width = 30, height = 25, dpi = 300, units = "cm")

########################## ILRUNgenesInf24hr ###################################################################### 
# 'Which' provides positions of 'TRUE' values.
sigILRUNgenesInf24hr<- ILRUNgenesInf24hr[ which( ILRUNgenesInf24hr$padj < alpha), ]
# Slices out rows where the adjusted p-value is <0.05.
sigILRUNgenesInf24hr <- sigILRUNgenesInf24hr[ which( abs(sigILRUNgenesInf24hr$log2FoldChange) > beta), ]
# Slices out rows where the fold change is above 2 and below -2.
# Abs=absolute, tells it to filter things above 2, ignoring whether value is positive or negative.
sigILRUNgenesInf24hr <- sigILRUNgenesInf24hr[ which(sigILRUNgenesInf24hr$baseMean > gamma), ]
write.csv(sigILRUNgenesInf24hr, file="results/sigILRUNgenesInf24hr_LT07_rev.csv")
# Slices out rows above an average count of 10 reads (anything below is rubbish).
ILRUNgenesInf24hr_hits <- rownames(sigILRUNgenesInf24hr)
length(ILRUNgenesInf24hr_hits)
write.csv(norm_counts[ILRUNgenesInf24hr_hits, ], file="results/DESeq2_sig_ILRUNgenesInf24hr_LT07_rev_normalized_counts.csv")
read_csv("results/sigILRUNgenesInf24hr_LT07_rev.csv") %>% 
  select("X1") %>% 
  rename(gene = X1) %>% 
  separate(gene, c("gene", "gene1")) %>% 
  select(gene) %>% 
  write.csv("results/ILRUNgenesInf24hr_LT07_rev_david.csv")

#making a PCA  
#read in the metadata
ILRUNgenesInf24hr_metadata <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv") %>% 
  filter(Infection == "Infected") %>% 
  filter(timepoint == "24hr") %>% 
  filter(Lane == "L001") %>%
  mutate(Name = str_extract(Sample, "^...."))

#join with the counts data 
ILRUNgenesInf24hr_expression <- read_csv("results/DESeq2_ilrun_siRNA_all_normalized_counts_LT07_rev.csv") %>% 
  rename(gene_locus = X1) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  right_join(ILRUNgenesInf24hr_metadata, by = "Sample") 

#scale gene expression for all samples 
ILRUNgenesInf24hr_expression_scaled_genes <- ILRUNgenesInf24hr_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -HTSeq_rev, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
ILRUNgenesInf24hr_pca_genes <- prcomp(ILRUNgenesInf24hr_expression_scaled_genes)

#tidy data frame for plotting
PCA_data_ILRUNgenesInf24hr_pca_genes <- ILRUNgenesInf24hr_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(ILRUNgenesInf24hr_metadata, by = "Sample") %>%
  spread(PC, expression)

plot_ILRUNgenesInf24hr <- ggplot(PCA_data_ILRUNgenesInf24hr_pca_genes, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = siRNA), size = 6)+
  geom_text(aes(label = Name), size = 4) +
  labs(title = "Principle Component Analysis ILRUNgenesInf24hr")

ggsave(filename = "results/PCA_plot_ILRUNgenesInf24hr.png", plot = plot_ILRUNgenesInf24hr, width = 20, height = 20, dpi = 300, units = "cm")
forplotting_ILRUNgenesInf6hr <- ILRUNgenesInf6hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

forplotting_ILRUNgenesInf24hr <- ILRUNgenesInf24hr %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>% 
  remove_rownames %>% 
  column_to_rownames(var="gene")

volcano_plot_ILRUNgenesInf24hr <- EnhancedVolcano(forplotting_ILRUNgenesInf24hr,
                                                 lab = rownames(forplotting_ILRUNgenesInf24hr),
                                                 x = 'log2FoldChange',
                                                 y = 'pvalue',
                                                 xlim = c(-3, 3),
                                                 title = 'ILRUN genes Infected 24hr',
                                                 pCutoff = 0.05,
                                                 FCcutoff = 0.75,
                                                 pointSize = 3.0,
                                                 labSize = 3.0)

ggsave(filename = "results/volcano_plot_ILRUNgenesInf24hr.png", plot = volcano_plot_ILRUNgenesInf24hr, width = 20, height = 20, dpi = 300, units = "cm")

############################################### Specific genes #######################################################################################
metadata_rev_LT07 <- read_csv("data/ilrun_metadata_siRNA_LT07_rev.csv")


#join with the counts data 
lessLT07_rev_expression <- read_csv("results/DESeq2_ilrun_siRNA_LT07_rev_all_normalized_counts.csv") %>%
  rename(gene_locus = X1) %>% 
  gather(Sample, expression, -gene_locus) %>%
  left_join(metadata_rev_LT07, by = "Sample") %>% 
  separate(gene_locus, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1) %>%
  mutate(Sample = str_remove(Sample, ""))

write.csv(lessLT07_rev_expression, file="results/DESeq2_ilrun_siRNA_LT07_rev_expression.csv")

tail(lessLT07_rev_expression)

TGFB <- read.csv("results/DESeq2_ilrun_siRNA_lessLT07_rev_expression.csv") %>%
  filter(gene_locus == "TGFB") %>% 
  ggplot(aes(x= Infection, y = expression, color = siRNA)) +
  geom_boxplot()+
  facet_wrap(~timepoint)+
  labs(title = "TGFB", 
       x = "") +
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 8 ),
         axis.title = element_text( size = 8))

ggsave(filename = "results/TGFB.png", plot = TGFB, width = 25, height = 20, dpi = 300, units = "cm")




##################################################ILRUN genes heat map##################################################################################

#tidy all the dataframes 
Uninf6hr <- read.csv("results/sigILRUNgenesUninf6hr_LT07_rev.csv") %>%
  as_tibble() %>% 
  rename(gene = X) %>% 
  separate(gene, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1, -baseMean, -lfcSE, -stat, -pvalue) %>% 
  rename(log2FC_Uninf6hr = log2FoldChange) %>% 
  rename(padj_Uninf6hr = padj)

  
Uninf24hr <- read.csv("results/sigILRUNgenesUninf24hr_LT07_rev.csv") %>%
  as_tibble() %>% 
  rename(gene = X) %>% 
  separate(gene, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1, -baseMean, -lfcSE, -stat, -pvalue) %>% 
  rename(log2FC_Uninf24hr = log2FoldChange) %>% 
  rename(padj_Uninf24hr = padj)


Inf6hr <- read.csv("results/sigILRUNgenesInf6hr_LT07_rev.csv") %>%
  as_tibble() %>% 
  rename(gene = X) %>% 
  separate(gene, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1, -baseMean, -lfcSE, -stat, -pvalue) %>% 
  rename(log2FC_Inf6hr = log2FoldChange) %>% 
  rename(padj_Inf6hr = padj)


Inf24hr <- read.csv("results/sigILRUNgenesInf24hr_LT07_rev.csv") %>%
  as_tibble() %>% 
  rename(gene = X) %>% 
  separate(gene, c("gene", "gene1"), sep = "_") %>% 
  select(-gene1, -baseMean, -lfcSE, -stat, -pvalue) %>% 
  rename(log2FC_Inf24hr = log2FoldChange) %>% 
  rename(padj_Inf24hr = padj)

#merge all the dataframes 

ILRUN_Uninf <- full_join(Uninf6hr, Uninf24hr, by = "gene" )

ILRUN_Inf <- full_join(Inf6hr, Inf24hr, by = "gene" )

ILRUN <-full_join(ILRUN_Uninf, ILRUN_Inf, by = "gene" )

#assigning the condition to a variable 

ILRUN_for_plotting <- ILRUN %>% 
  select(gene, log2FC_Uninf6hr, log2FC_Uninf24hr, log2FC_Inf6hr, log2FC_Inf24hr) %>%
  mutate(FC_6hf = log2FC_Inf6hr - log2FC_Uninf6hr) %>% 
  mutate(FC_24hr = log2FC_Inf24hr - log2FC_Uninf24hr) %>% 
  filter(FC_6hf >= 0.25 | FC_6hf <= -0.25 | FC_24hr >=0.25 | FC_24hr < -0.25 ) %>%
  filter(str_detect(gene, "LOC") == FALSE) %>% 
  filter(FC_24hr != "NA") %>% 
  select(-FC_6hf) %>% 
  select(-FC_24hr) %>% 
  gather(condition, log2FC,  -gene) %>% 
  separate(condition, c("log", "condition"), sep = "_") %>% 
  select(-log) %>% 
  mutate(direction = log2FC>0) %>% 
  filter(log2FC != "NA")
  

# Negative numbers are increasingly red, positive numbers are increasingly green
# Transitions through white at 0 (has a 'midpoint' argument but default is 0)
# Can use colour names found in colours() or can use hexcodes eg "#ff0000" = red
scale_fill_gradient2(low = "red", mid = "white", high = "green")


ILRUN_heatmap <- ggplot(ILRUN_for_plotting, aes(x = condition, y=gene, fill=log2FC)) + 
  geom_tile()+
  scale_fill_distiller(palette = "PRGn", limits = c(-3, 3))
  
  


# Converts one of the colorbrewer2.org palettes for use with continuous data
# The diverging palettes from there are good for heatmaps
# Eg, for the purple/green palette find the name from the site (PRGn)
# Doesn't have a midpoint setting, so you'll need to set the limits manually to 
# be symmetric based on the largest absolute log2FC
scale_fill_distiller(palette = "PRGn", limits = c(-3, 3)
                     
                     
  
  theme(axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 10))+
  labs(title = "ILRUN DE genes altered by Imfection", color = "ILRUN upregulation") +
  labs
                            

ggsave(filename = "results/ILRUN_heatmap.png", plot = ILRUN_heatmap, width = 10, height = 30, dpi = 300, units = "cm")
