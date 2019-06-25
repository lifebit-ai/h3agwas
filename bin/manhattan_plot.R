#!/usr/bin/env Rscript

library(manhattanly)
library(plotly)

args = commandArgs(trailingOnly=TRUE)

assoc_file <- args[1] # eg "~/Documents/GitHub/lifebit-gwas/manhattan-plot/header.manhattan_plot_data.assoc"
assoc <- read.table(assoc_file, quote="\"", comment.char="", header=T)
# plot <- manhattanly(assoc, snp="SNP", chr="CHR")
# plot <- manhattanly(subset(assoc, CHR %in% 20:21), snp="SNP", chr="CHR")
# htmlwidgets::saveWidget(as.widget(plot), "manhattan_plot.html")

small_assoc <- assoc[assoc$CHR == 21,]
plot <- manhattanly(small_assoc, snp="SNP", chr="CHR")
htmlwidgets::saveWidget(as.widget(plot), "manhattan_plot.html")