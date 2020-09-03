##############################################################################
###                                                                        ###
###              Reads mapping to viral genome Coverview                   ###
###                                                                        ###
##############################################################################
#coverview_virus.r

# remind R where to look for libraries
.libPaths(c("C:/Users/ale097/Data School/Packages"))
# load libraries
library(tidyverse)
library(dplyr)
library(readr)
library(vroom)
library(stringr)


###First get the list of all txt files
list_of_coverview_summary_files <- list.files(path = ".", recursive = TRUE,
                            pattern = "\\.txt$", 
                            full.names = TRUE)

####total read counts for normalisation
total_counts <- read.csv("results/total_counts.csv") %>% 
  as_tibble() %>%
  filter(Sample != "LT41") %>% 
  filter(Sample != "LT42") %>%
  filter(Sample != "LT21") %>% 
  filter(Sample != "LT22") %>% 
  filter(Sample != "LT23") %>%
  filter(Sample != "LT24") %>%
  filter(Sample != "LT22") %>% 
  filter(Sample != "LT23") %>%
  filter(Sample != "LT24") %>%
  select(Reads)

#read in the rRNA data
counts_18S <- read.csv("results/counts_18s.csv")

#read in the metadata
metadata_viral_readcount <- read.csv("data/ilrun_metadata_siRNA.csv") %>% 
  rename(Sample_Lane = Sample) %>% 
  filter(Sample_Lane !="LT21_L001") %>%
  filter(Sample_Lane !="LT21_L002") %>%
  filter(Sample_Lane !="LT22_L001") %>%
  filter(Sample_Lane !="LT22_L002") %>%
  filter(Sample_Lane !="LT23_L001") %>%
  filter(Sample_Lane !="LT23_L002") %>%
  filter(Sample_Lane !="LT24_L001") %>%
  filter(Sample_Lane !="LT24_L002") 


###Put them into a dataframe using vroom
coverview_ilrun_SARSCOV2 <- vroom(list_of_coverview_summary_files, id = "FileName") %>% 
  filter(`#CHROM` == "MT007544.1" ) %>% 
  mutate(Sample_Lane = str_extract(FileName, "LT.......")) %>% 
  select(-RCIN, -RCOUT, -FileName ) %>% 
  mutate(Sample = str_extract(Sample_Lane, "LT..")) %>% 
  mutate(Lane = str_extract(Sample_Lane, "L00.")) %>%
  left_join(metadata_viral_readcount, by = "Sample_Lane") %>% 
  select(-Sample_Lane, -Lane.y) %>% 
  filter(Infection != "NA") %>% 
  bind_cols(total_counts) %>% 
  mutate(cpm = RC/Reads*1000000 )

viral_reads_ilrun <- ggplot(coverview_ilrun_SARSCOV2, aes(x = Infection, y = cpm, color = siRNA)) +
  geom_boxplot()+
  scale_y_log10()+
  labs(title = "SARS-CoV2 read counts")

ggsave(filename = "results/viral_reads_ilrun_july.png", plot = viral_reads_ilrun, width = 12, height = 10, dpi = 300, units = "cm")

#want to just plot the infected samples

coverview_ilrun_SARSCOV2_infected <- coverview_ilrun_SARSCOV2 %>% 
  filter(Infection == "Infected " )

write.csv(coverview_ilrun_SARSCOV2_infected, "results/coverview_ilrun_SARSCOV2_infected_cpm.csv")

viral_reads_ilrun_infected <- ggplot(coverview_ilrun_SARSCOV2_infected, aes(x = timepoint, y = cpm, color = siRNA)) +
  geom_boxplot()+
  scale_y_log10()+
  labs(title = "SARS-CoV2 read counts for infected samples")


ggsave(filename = "results/viral_reads_ilrun_infected.png", plot = viral_reads_ilrun_infected, width = 12, height = 10, dpi = 300, units = "cm")



#######################18S rRNA contamination######################

total_counts_for18S <- read.csv("results/total_counts.csv") %>% 
  as_tibble() %>%
  filter(Sample != "LT41") %>% 
  filter(Sample != "LT42")

counts_18S <- read.csv("results/counts_18s.csv") %>% 
  as_tibble() %>%
  filter( X.CHROM != "LT41_L001") %>% 
  filter( X.CHROM != "LT41_L002") %>% 
  filter( X.CHROM != "LT42_L001") %>% 
  filter( X.CHROM != "LT42_L002") %>%  
  bind_cols(total_counts_for18S) %>% 
  mutate(cpm = RC/Reads*100000) %>% 
  rename(Sample_Lane = X.CHROM) %>%
  left_join(metadata_viral_readcount, by = "Sample_Lane") %>% 
  filter(Lane == "L001") %>% 
  mutate(Sample = str_remove(Sample, "LT"))

  
  
rRNA_contamination <- ggplot(counts_18S, aes(x = Sample, y = cpm, color = Infection )) +
  geom_point(size = 4) +
  labs(title = "18S rRNA contamination")

ggsave(filename = "results/rRNA_contamination.png", plot = rRNA_contamination, width = 22, height = 10, dpi = 300, units = "cm")


