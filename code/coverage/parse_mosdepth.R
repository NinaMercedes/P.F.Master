library(dplyr)
library(ggplot2)
library(showtext)
library(readr)
showtext_auto()

# wd <- "/mnt/storage9/emilia/Psimium/analysis/nanopore_assembly/map_PvP01/minimap2/coverage"
# samples <- read.table("/mnt/storage9/emilia/Psimium/data/nanopore/PvPs_nanopore_samples.txt")$V1

wd <- "."
samples <- read.table("samples.txt")$V1

summary <- c()
i <- 0
for (sample in samples) {
    cov <- read_tsv(file.path(wd, sprintf("%s.regions.bed.gz", sample)), col_names = FALSE) %>%
        mutate(X2 = X2 / 500, X3 = X3 / 500)
    sum <- read_tsv(file.path(wd, sprintf("%s.mosdepth.summary.txt", sample)))
    sumdf <- sum %>%
     select(chrom, mean) %>%
     filter(!grepl("_region", chrom)) %>% rename(!!sym(sample) := "mean")

    if (i == 0) {
        summary <- sumdf
    } else {
        summary <- summary %>% inner_join(sumdf)
    }
    pdf(file.path(wd, paste0(sample, "_cov_per_chr.pdf")), width = 10, height = 4)
    for (chr in unique(cov$X1)[-c(15, 16)]) {
        cov_sub <- cov %>% filter(X1 == chr)
        sum_sub <- sum %>% filter(chrom == chr)
        pl <- ggplot(cov_sub, aes(x = X3, y = X4)) +
            geom_line() +
            theme_bw() +
            labs(x = "Genomic window (1kb)",
                y = "Mean Coverage",
                title = chr,
                subtitle = sprintf("mean: %.2f | min: %d | max: %d", sum_sub$mean, sum_sub$min, sum_sub$max))
        print(pl)
    }
    i = i + 1
    dev.off()
}

write.table(summary, file.path(wd, "summary_coverage_new_batch.txt"), quote=F, row.names=F, col.names=T, sep="\t")
