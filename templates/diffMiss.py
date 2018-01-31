#!/usr/bin/env python3

#Load SNP differential missingness file and generate distribution
args <- commandArgs(TRUE)
b.frq <- read.table("$input",header=T)
if (nrow(b.frq) >= 1) {
b.frq\$logP = log10(b.frq\$P)
pdf("$output")
plot(ecdf(b.frq\$logP), xlim=c(-10,0),ylim=c(0,1),pch=20, main="Distribution of differential missingness P-values", xlab="logP Differential Missingness", ylab="Fraction of SNPs",col="red",axes=T)
} else {
    print("No differential missingness info to plot")}
