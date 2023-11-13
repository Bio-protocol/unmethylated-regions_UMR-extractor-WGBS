#!/usr/bin/Rscript
##########

#Peter Crisp
#Leroy Mangila
#2022-06-17
#R script to classify methylation levels in 100bp tile methylation data

# Notes
# We suggest removing organelles from the data before proceeding with this step

# Arguments
## Defaults
input   <- "./"
out_dir <- "./"

#Sample selection #####
args <- commandArgs(trailingOnly=TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned
reference_100bp_tiles <- args[1]
coverage_filter_min   <- as.double(args[2])
site_filter_min       <- as.double(args[3])
MR_percent            <- as.double(args[4])
UMR_percent           <- as.double(args[5])
input                 <- as.character(args[6])
out_dir               <- as.character(args[7])

######## argument recommendations
######## de bug
# args
#reference_100bp_tiles = "maize_chr1_reference_100bp_tiles_sites_counts.txt"
#coverage_filter_min = 3
#site_filter_min = 2
#MR_percent = 0.4
#UMR_percent = 0.1

# The following arguments are required in this order, in bracket are the arguments for this example:
# The reference genome site file (“maize_chr1_reference_100bp_tiles_sites_counts.txt”)
# Minimum coverage (suggestion 3x or 5x "3")
# Minimum number of sites (suggestion 2 "2")
# Minimum percent to be considered methylated (suggestion 40% "0.4")
# Maximum percent to be considered unmethylated (suggestion 10% "0.1")


###########################
# load required libraries
if (!require(tidyverse)) {install.packages("tidyverse")}
library(tidyverse)
old.scipen <- getOption("scipen")
options(scipen=999)

text_size_theme_8 <- theme(axis.text=element_text(size=8),
                           axis.title=element_text(size=8),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.title=element_text(size=8),
                           legend.text=element_text(size=8))

###########################
# Module #1
###########################

#############
## read in mC tile data
#############

reference_tiles <- read_tsv(reference_100bp_tiles, col_names = T,
                            cols(
                              chr = col_character(),
                              start = col_integer(),
                              end = col_integer()))

# Print total size of the reference genome
reference_tiles %>% mutate(size = end - start) %>% summarise(MB = sum(size)/1000000)

###########
## CG  ##
# read in CG and parse
bedGraph <- read_tsv(paste0(input, "/BSMAP_out.txt.100.CG.bed"), col_names = F,
                            cols(
                              X1 = col_character(),
                              X2 = col_integer(),
                              X3 = col_integer(),
                              X4 = col_number()
                            ))

colnames(bedGraph) <- c("chr", "start_zBased", "end", "sites_with_data", "C", "CT", "ratio")

# remove NAs - this is necessary because the output from the perl script
# is any tile with data in any context per sample
bedGraph <- na.omit(bedGraph)

#convert to one-based coordinate and sort
CG <-bedGraph %>% 
  mutate(start = start_zBased +1) %>%
  left_join(reference_tiles, by = c("chr", "start", "end")) %>%
  select(chr, start, end, C, CT, ratio, sites_with_data, cg_sites) %>% arrange(chr, start) %>% 
  mutate(cov = CT/cg_sites)


###########
## CHG  ##
# read in CHG and parse
bedGraph <- read_tsv(paste0(input, "/BSMAP_out.txt.100.CHG.bed"), col_names = F,
                     cols(
                       X1 = col_character(),
                       X2 = col_integer(),
                       X3 = col_integer(),
                       X4 = col_number()
                     ))

colnames(bedGraph) <- c("chr", "start_zBased", "end", "sites_with_data", "C", "CT", "ratio")
# Remove NAs - this is necessary because the output from the perl script
# is any tile with data in any context per sample
bedGraph <- na.omit(bedGraph)

# convert to one-based coordinate and sort
CHG <-bedGraph %>% 
  mutate(start = start_zBased +1) %>%
  left_join(reference_tiles, by = c("chr", "start", "end")) %>%
  select(chr, start, end, C, CT, ratio, sites_with_data, chg_sites) %>% arrange(chr, start) %>% 
  mutate(cov = CT/chg_sites)

###########
## CHH  ##
# read in CHH and parse
bedGraph <- read_tsv(paste0(input,"/BSMAP_out.txt.100.CHH.bed"), col_names = F,
                     cols(
                       X1 = col_character(),
                       X2 = col_integer(),
                       X3 = col_integer(),
                       X4 = col_number()
                     ))

colnames(bedGraph) <- c("chr", "start_zBased", "end", "sites_with_data", "C", "CT", "ratio")
# Remove NAs - this is necessary because the output from the perl script
# is any tile with data in any context per sample
bedGraph <- na.omit(bedGraph)

#convert to one-based coordinate and sort
CHH <-bedGraph %>% 
  mutate(start = start_zBased +1) %>%
  left_join(reference_tiles, by = c("chr", "start", "end")) %>%
  select(chr, start, end, C, CT, ratio, sites_with_data, chh_sites) %>% arrange(chr, start) %>% 
  mutate(cov = CT/chh_sites)

########### ########### ###########
## Filter
########### ########### ###########

# filter
CG_ratio <- CG %>%
  mutate(CG = ifelse(cg_sites < site_filter_min | cov < coverage_filter_min, NA, ratio)) %>%
  select(chr, start, end, CG)

# filter
CHG_ratio <- CHG %>%
  mutate(CHG = ifelse(chg_sites < site_filter_min | cov < coverage_filter_min, NA, ratio)) %>%
  select(chr, start, end, CHG)

# filter
CHH_ratio <- CHH %>%
  mutate(CHH = ifelse(chh_sites < site_filter_min | cov < coverage_filter_min, NA, ratio)) %>%
  select(chr, start, end, CHH)

############
# merge with reference
# call tiles with no sites
merged_mC <- reference_tiles %>%
  # select("chr", "start", "end") %>%
  left_join(CG_ratio, by = c("chr", "start", "end")) %>%
  left_join(CHG_ratio, by = c("chr", "start", "end")) %>%
  left_join(CHH_ratio, by = c("chr", "start", "end")) %>%
  mutate(cg_sites = ifelse(cg_sites < site_filter_min, "n", "y"),
         chg_sites = ifelse(chg_sites < site_filter_min, "n", "y"),
         chh_sites = ifelse(chh_sites < site_filter_min, "n", "y"))

# print the top of each file to check the filtered files were created
print(reference_tiles, n = 20)
print(CHG_ratio, n = 20)
print(merged_mC, n = 20)


# strongly recommended to remove organelles, below code is commented out as it is specific for maize
# merged_mC %>% distinct(chr)
# merged_mC %>% distinct(chr) %>% filter(grepl(paste("M", "P", sep = "|"), chr))
# merged_mC_sans_orgs <- merged_mC %>% filter(!chr %in% c("Mt", "Pt"))

merged_mC_sans_orgs <- merged_mC

###################
# distribution and averages of mC levels
# uncomment to run these summary plots

# # distro
# plot_data <- merged_mC_sans_orgs %>% slice(1:10000) %>%
#   select(CG:CHH) %>%
#   gather(key = context, value = percent) %>%
#   mutate(percent = percent *100)
# 
# g <- ggplot(plot_data, aes(x = percent)) +
#   geom_density() +
#   geom_vline(xintercept = c(UMR_percent*100, MR_percent*100), colour = 'blue', linetype = 'dashed') +
#   facet_grid(context ~., scales = 'free') +
#   theme_minimal() +
#   text_size_theme_8
# 
# ggsave(plot = g, filename = paste0(out_dir, "/mC_tile_density.pdf"), h = 4, w = 2)
# 
# # mean
# plot_data_average <- plot_data %>%
#   group_by(context) %>%
#   summarise(mean = mean(percent, na.rm = T))
# 
# g <- ggplot(plot_data_average, aes(x = context, y = mean)) +
#   geom_bar(stat = 'identity') +
#   ylim(0,100) +
#   theme_minimal() +
#   text_size_theme_8
# 
# ggsave(plot = g, filename = paste0(out_dir, "/mC_tile_average_mC.pdf"), h = 2, w = 1.1)

# write.table(plot_data_average, paste0(out_dir,"/", sample_to_crunch, "mC_tile_average_mC.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

###################
## Annotate domains

### CHG mC domains
mC_domains <- merged_mC_sans_orgs %>%
  mutate(domain_tmp = ifelse(chg_sites == "n", "no_sites", #this catches tiles with no CHG sites
                             ifelse(CHG >= MR_percent, "MR", #this catches sites with no data too: they get NA
                                    ifelse(CHG < UMR_percent, "UMR", "Intermediate")))) %>%
  mutate(CHG_based_domain = ifelse(is.na(domain_tmp), "Missing_Data", domain_tmp))

mC_domains


### 40% mC
# be careful this steps are maize specific and would need to be adjusted to another species...
mC_domains2 <- mC_domains %>%
  mutate(domain_tmp = ifelse(chh_sites == "n", "no_sites", #this will catch the 1.45% of the genome that lack CHH sites (or any other site)
                             ifelse(CHH >= 0.15, "RdDM", #
                                    ifelse(cg_sites == "n" & chg_sites == "y" & CHG >= MR_percent, "Heterochromatin", # this catches tiles with CHG mC but no CG sites (otherwise they would get NA)
                                           ifelse(cg_sites == "y" & chg_sites == "y" & CHG >= MR_percent & CG >= MR_percent, "Heterochromatin",
                                                  ifelse(cg_sites == "y" & CG >= MR_percent, "CG_only",
                                                         ifelse(cg_sites == "y" & chg_sites == "y" & CG < UMR_percent & CHG < UMR_percent & CHH < UMR_percent, "Unmethylated",
                                                                ifelse(cg_sites == "n" & chg_sites == "y" & CHG < UMR_percent & CHH < UMR_percent, "Unmethylated",
                                                                       ifelse(cg_sites == "y" & chg_sites == "n" & CG < UMR_percent & CHH < UMR_percent, "Unmethylated",
                                                                              ifelse(cg_sites == "y" & chg_sites == "y" & CG >= UMR_percent & CHG >= UMR_percent & CHH >= UMR_percent, "Intermediate",
                                                                                     ifelse(cg_sites == "n" & chg_sites == "y" & CHG >= UMR_percent & CHH >= UMR_percent, "Intermediate",
                                                                                            ifelse(cg_sites == "y" & chg_sites == "n" & CG >= UMR_percent & CHH >= UMR_percent, "Intermediate",
                                                                                                   ifelse(cg_sites == "n" | chg_sites == "n", "no_sites", NA))))))))))))) %>%
  mutate(domain = ifelse(is.na(domain_tmp), "Missing_Data", domain_tmp)) %>%
  mutate(domain_simple = ifelse(domain == "Heterochromatin", "MR",
                                ifelse(domain == "Missing_Data", "no_data",
                                       ifelse(domain == "Unmethylated", "UMR",
                                              ifelse(domain == "no_sites", "No_sites","other_mC")))))

mC_domains2

############
## summarise

mC_domains_freq <- mC_domains2 %>%
  group_by(domain) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  mutate(percent = total/sum(total)*100,
         MB = total*100/1000000)
mC_domains_freq

write.table(mC_domains_freq, paste0(out_dir, "/mC_domains_freq.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

######################### #########################

######################### #########################
# write bed file
mC_domains_bed <- mC_domains2 %>%
  mutate(start = start-1,
         score = ".",
         strand = ".") %>%
  select(chr, start, end, domain, score, strand)


# write the whole data file
write.table(mC_domains2, paste0(out_dir, "/mC_domains",
                                "_cov_",coverage_filter_min,
                                "_sites_",site_filter_min,
                                "_MR_",MR_percent,
                                "_UMR_",UMR_percent,".txt"), sep = "\t", quote = F, row.names = F, col.names = T)


## Make UMT and ND only bed files
# make UMR only bedfile
mC_domains2 %>% distinct(chr)

############# UMTs
# subset to UMTs
mC_domains2 %>% distinct(domain)
UMT_only <- mC_domains2 %>%
  filter(domain == "Unmethylated") %>%
  mutate(start = start-1) %>%
  select(chr:end, domain)
UMT_only
# 1,071,742

write.table(UMT_only, paste0(out_dir, "/UMTs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

############# Missing data and no sites
# subset to NDs
ND_only <- mC_domains2 %>%
  filter(domain %in% c("Missing_Data", "no_sites")) %>%
  mutate(start = start-1) %>%
  select(chr:end, domain)
ND_only
# 3,372,705

write.table(ND_only, paste0(out_dir, "/NDs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

############# Tiles with data (exclude no data and no sites)

NDs_data <- mC_domains_bed %>%
filter(!domain %in% c("no_sites", "Missing_Data")) %>%
select(chr, start, end, domain)

write.table(NDs_data,
            paste0(out_dir, "/mC_domains",
                   "_cov_",coverage_filter_min,
                   "_sites_",site_filter_min,
                   "_MR_",MR_percent,
                   "_UMR_",UMR_percent, "_tiles_with_data",".bed"),
            sep = "\t", quote = F, row.names = F, col.names = F)
