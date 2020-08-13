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
write.csv(norm_counts, file="results/DESeq2_ilrun_siRNA_all_normalized_counts.csv")


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



#####Repeat analysis for infection versus uninfected counts dataset excluding LT07 and excluding lane duplicates for plot readability###############

# Construct a SummarizedExperiment object:
dds_inf <- DESeqDataSetFromMatrix(
  countData = df,
  colData = md,
  design = ~ Infection + Lane) # ~ is representative of 'by', i.e. compare by condition, sex.

# Perform DE testing:
dds_inf <- DESeq(dds_inf)

# Output normalized counts:
norm_counts_inf <- counts (dds_inf, normalized=TRUE)
write.csv(norm_counts_inf, file="results/DESeq2_ilrun_siRNA_inf_normalized_counts.csv")

#read in the metadata
Inf_lessLT07_lessLane_metadata <- read_csv("data/ilrun_metadata_siRNA.csv") %>% 
  mutate(Name = str_extract(Sample, "^....")) %>%
  filter(Sample != "LT07_L001" |Sample != "LT07_L002") %>%
  filter(Lane != "L002")

#join with the counts data 
Inf_lessLT07_lessLane_expression <- read_csv("results/DESeq2_ilrun_siRNA_inf_normalized_counts.csv") %>%
  rename(gene_locus = X1) %>%
  select(-LT07_L001, -LT07_L002) %>%
  gather(Sample, expression, -gene_locus) %>%
  filter(str_detect(Sample, 'L001')) %>% 
  left_join(lessLT07_lessLane_metadata, by = "Sample")

#scale gene expression for all samples 
Inf_lessLT07_lessLane_scaled_genes <- Inf_lessLT07_lessLane_expression %>%
  spread(gene_locus, expression) %>%
  select(-Infection, -siRNA, -timepoint, -Condition, -Lane, -Name) %>% 
  column_to_rownames("Sample") %>% 
  scale()

#use the prcomp function on scaled samples
Inf_lessLT07_Lesslane_pca_genes <- prcomp(Inf_lessLT07_lessLane_scaled_genes)

#tidy data frame for plotting
PCA_data_lessLT07_Lesslane_inf <- Inf_lessLT07_Lesslane_pca_genes$x %>% 
  as_tibble(rownames = "Sample") %>%
  gather(PC, expression, -Sample) %>% 
  left_join(all_metadata, by = "Sample") %>%
  spread(PC, expression)

text_plot_lessLT07_Lesslane_inf <- ggplot(PCA_data_lessLT07_Lesslane_inf, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = siRNA, color = Infection), size = 6)+
  geom_text(aes(label = timepoint), size = 4) +
  labs(title = "Principle Component Analysis less LT07 less Lanes infection normalized counts")

ggsave(filename = "results/PCA_lessLT07_LessLane_inf.png", plot = text_plot_lessLT07_Lesslane_inf, width = 15, height = 10, dpi = 300, units = "cm")


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



PDZD4 <- lessLT07_expression %>%
  filter(gene_locus == "XLOC_037006_PDZD4") %>% 
  ggplot(aes(x= siRNA, y = expression)) +
  geom_boxplot()+
  labs(title = "PDZD4", 
       x = "")
