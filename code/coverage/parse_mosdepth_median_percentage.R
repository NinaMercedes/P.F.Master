library(dplyr)
library(ggplot2)
#library(showtext)
library(readr)
library(data.table)
#showtext_auto()

# wd <- "/mnt/storage9/emilia/Psimium/analysis/nanopore_assembly/map_PvP01/minimap2/coverage"
# samples <- read.table("/mnt/storage9/emilia/Psimium/data/nanopore/PvPs_nanopore_samples.txt")$V1

wd <- "."
samples <- read.table("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/new_vcf.txt")$V1
threshold1 <- 5
threshold2 <- 10
threshold3 <- 15
threshold4 <- 20
threshold5 <- 30
percent <- list()
cov <-list()
for (i in 1:length(samples)){
    if (file.exists(file.path(wd, sprintf("%s.regions.bed.gz", samples[i])))){
        cas <- read_tsv(file.path(wd, sprintf("%s.regions.bed.gz", samples[i])), col_names = FALSE) 
        if (nrow(cas)>0){
            cov[[i]] <- cas
            colnames(cov[[i]]) <- c("Chromosome", "Start", "End", "Coverage")
            median_cov <- median(cov[[i]]$Coverage)
            min_cov <- min(cov[[i]]$Coverage)
            max_cov <- max(cov[[i]]$Coverage)
            greater_than1 <- cov[[i]] %>% filter(Coverage>threshold1)
            greater_than2 <- cov[[i]] %>% filter(Coverage>threshold2)
            greater_than3 <- cov[[i]] %>% filter(Coverage>threshold3)
            greater_than4 <- cov[[i]] %>% filter(Coverage>threshold4)
            greater_than5 <- cov[[i]] %>% filter(Coverage>threshold5)
            cov[[i]]$sample_id <- samples[[i]]
            percent[[i]]<- data.frame(samples[[i]], (nrow(greater_than1)*100)/nrow(cov[[i]]),(nrow(greater_than2)*100)/nrow(cov[[i]]),(nrow(greater_than3)*100)/nrow(cov[[i]]),(nrow(greater_than4)*100)/nrow(cov[[i]]),(nrow(greater_than5)*100)/nrow(cov[[i]]), median_cov, min_cov, max_cov)
            colnames(percent[[i]]) <- c("sample", "X5 percentage", "X10 percentage","X15 percentage","X20 percentage","X30 percentage", "Median Coverage", "Minimum Coverage", "Maximum Coverage")
        }
    }
}
percentages <- rbindlist(percent)
percentages$PASS <- ifelse(percentages$`X5 percentage`>60, TRUE, FALSE)
#coverages <- rbindlist(cov)
#write.csv(coverages, "coverage_per_region.csv", row.names=FALSE)
write.csv(percentages, "/mnt/storage13/nbillows/Pf_09_24/coverage_new_samp.csv", row.names=FALSE)


