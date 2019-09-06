#!/usr/bin/env Rscript

library(qqman)
# library(manhattanly)
# library(plotly)

args = commandArgs(trailingOnly=TRUE)
assoc_file <- args[1] # eg "out.assoc"
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line
blues.c <- c("#6E65C2", "#4D43AE", "#3226A6", "#211785", "#150D69")

assoc <- read.table(assoc_file, quote="\"", comment.char="", header=T)
assoc <- assoc[!is.na(assoc$P),]

# subset
# small_assoc <- assoc[assoc$CHR == 21,]
# interactive HTML plot
# plot <- manhattanly(subset(assoc, CHR %in% 20:21), snp="SNP", chr="CHR")
# htmlwidgets::saveWidget(as.widget(plot), "manhattan_plot.html")

workdir <- getwd()
png_name <- paste0(workdir, "/manhattan_plot.png")
png(png_name, width=1425, height=975)
print(manhattan(assoc, 
                suggestiveline = -log10(sugg),
                genomewideline = -log10(sig),
                col=blues.c, 
                cex=1.1,
                cex.axis = 2,
                cex.lab = 1.5,
                cex.main = 2,
                annotatePval = sugg, annotateTop = FALSE,
                main = "Manhattan Plot GWAS Results\nQQman Plot"))
dev.off()

