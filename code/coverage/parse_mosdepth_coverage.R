library(dplyr)
library(ggplot2)
library(showtext)
library(readr)
library(data.table)
showtext_auto()

# wd <- "/mnt/storage9/emilia/Psimium/analysis/nanopore_assembly/map_PvP01/minimap2/coverage"
# samples <- read.table("/mnt/storage9/emilia/Psimium/data/nanopore/PvPs_nanopore_samples.txt")$V1

wd <- "."
samples <- read.table("samples.txt")$V1
threshold <- 5

percent <- list()
for (i in 1:length(samples)) {
    cov <- read_tsv(file.path(wd, sprintf("%s.regions.bed.gz", samples[i])), col_names = FALSE) %>%
        mutate(X2 = X2 / 500, X3 = X3 / 500)
    greater_than <- cov %>% filter(X4>threshold)
    percent[[i]]<- data.frame(samples[[i]], (nrow(greater_than)*100)/nrow(cov))
    colnames(percent[[i]]) <- c("sample", "X5 percentage")
}
percentages <- rbindlist(percent)
percentages$PASS <- ifelse(percentages$`X5 percentage`>60, TRUE, FALSE)
write.csv(percentages, "coverage_percentage.csv", row.names=FALSE)


