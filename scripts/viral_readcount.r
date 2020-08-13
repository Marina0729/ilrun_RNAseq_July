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


###First get the list of all txt files
list_of_coverview_summary_files <- list.files(path = ".", recursive = TRUE,
                            pattern = "\\.txt$", 
                            full.names = TRUE)
#read in the metadata
metadata_viral_readcount <- read.csv("data/ilrun_metadata_siRNA.csv") %>% 
  rename(Sample_Lane = Sample)


###Put them into a dataframe using vroom
coverview_ilrun_SARSCOV2 <- vroom(list_of_coverview_summary_files, id = "FileName") %>% 
  filter(`#CHROM` == "MT007544.1" ) %>% 
  mutate(Sample_Lane = str_extract(FileName, "LT.......")) %>% 
  select(-RCIN, -RCOUT, -FileName ) %>% 
  mutate(Sample = str_extract(Sample_Lane, "LT..")) %>% 
  mutate(Lane = str_extract(Sample_Lane, "L00.")) %>%
  left_join(metadata_viral_readcount, by = "Sample_Lane") %>% 
  select(-Sample_Lane, -Lane.y) %>% 
  filter(Infection != "NA")

viral_reads_ilrun <- ggplot(coverview_ilrun_SARSCOV2, aes(x = Infection, y = RC, color = siRNA)) +
  geom_boxplot()+
  scale_y_log10()+
  labs(title = "SARS-CoV2 read counts")

ggsave(filename = "results/viral_reads_ilrun_july.png", plot = viral_reads_ilrun, width = 12, height = 10, dpi = 300, units = "cm")

#want to just plot the infected samples

coverview_ilrun_SARSCOV2_infected <- coverview_ilrun_SARSCOV2 %>% 
  filter(RC > 1000 )

viral_reads_ilrun_infected <- ggplot(coverview_ilrun_SARSCOV2_infected, aes(x = timepoint, y = RC, color = siRNA)) +
  geom_boxplot()+
  scale_y_log10()+
  labs(title = "SARS-CoV2 read counts")


ggsave(filename = "results/viral_reads_ilrun_infected.png", plot = viral_reads_ilrun_infected, width = 12, height = 10, dpi = 300, units = "cm")
