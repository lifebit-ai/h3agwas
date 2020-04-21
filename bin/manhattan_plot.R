#!/usr/bin/env Rscript

# install.packages("CMplot")
library(qqman)
library(plotly)
library(tidyverse)
library(CMplot)
library(manhattanly)

args = commandArgs(trailingOnly=TRUE)
assoc_file <- args[1] # eg "out.assoc"
annot_file <- args[2] # eg "annotation_GWAS_subset.csv"

##########
### qqman
##########
assoc <- read.table(assoc_file, quote="\"", comment.char="", header=T)
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line
blues.c <- c("#6E65C2", "#4D43AE", "#3226A6", "#211785", "#150D69")
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

##############################
### Plotly for Manhattan + QQ 
##############################

annotation <- read.csv(annot_file, stringsAsFactors = F)
geno <- annotation %>% select(c(rsid, Gene, Function))  %>%
  mutate_at(vars(Function), list(~ifelse(.=="NA", "intergenic",.)))
gwas_input <- read.table(assoc_file, h=T, stringsAsFactors = F) %>% left_join(., geno, by=c("SNP"="rsid")) %>% mutate_at(vars(CHR,BP, OR),list(~replace_na(.,"")))%>% mutate_at(vars(CHR,BP, P, OR), list(~as.numeric(.))) 

# Prepare the dataset
gwas_plot <- gwas_input %>% 
  group_by(CHR) %>%    # Compute chromosome size
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
  select(-chr_len) %>%
  left_join(gwas_input, ., by=c("CHR"="CHR")) %>% # Add this info to the initial dataset
  arrange(CHR, BP) %>% # Add a cumulative position of each SNP
  mutate(BPcum=BP+tot) %>%
  mutate_at(vars(Gene,Function),list(~replace_na(.,""))) %>%
  mutate(volcano_threshold=case_when(
    .$OR <= -3 | .$OR >= 3 ~ "0",
    TRUE ~ "1",
  ))

# Prepare X axis
axisdf <- gwas_plot %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Prepare text description for each SNP:
gwas_plot$text <- paste("SNP: ", gwas_plot$SNP, "\nPosition: ", gwas_plot$BP, "\nChromosome: ", gwas_plot$CHR, "\nREF: ", gwas_plot$A1, "\nALT: ", gwas_plot$A2,
                        "\nFunction: ",gwas_plot$Function,"\nMAF: ", gwas_plot$F_A, sprintf("\nP-value: %2e", round(gwas_plot$P,3)),
                        "\nGene: ",gwas_plot$Gene, "\nOdds Ratio: ", gwas_plot$OR, sep="")

# Manhattan
p <- ggplot(gwas_plot, aes(x=BPcum, y=-log10(P), text=text)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  ggtitle("Manhattan Plot") +
  xlab("Chromosome") + 
  ylab("-log10 p-value") +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

manhattan <- ggplotly(p, tooltip="text")

# Volcano
p1 <-ggplot(gwas_plot, aes(x=OR, y=-log10(P), text=text)) +
  geom_point(aes(colour=volcano_threshold), alpha=0.8, size=1.3) +
  xlab("Odds Ratio") + 
  ylab("-log10 p-value") +
  ggtitle("Volcano Plot") +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "skyblue", size=1)+
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "skyblue", size=0.5)+
  scale_color_manual(values = rep(c("red", "skyblue"), 22 )) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
  )

volcano <- ggplotly(p1, tooltip="text")

qqplot <- qqly(gwas_plot %>% filter(P > 0), main= "QQ plot", snp="SNP", gene = "Gene")

#qqly(HapMap, snp = "SNP", gene = "GENE")

# Annotation table
top_markers <- gwas_plot %>% dplyr::arrange(P) %>% slice(1:500) %>% select(SNP)
annotation_top <- annotation %>% select(
  Location,rsid,Reference,Alternative,Gene,Function,SIFT.Effect,Polyphen2_HDIV.Effect,MutationTaster.Effect, CADD.Score,
  gnomAD.AF.All, gnomAD.AF.NFE) %>% filter(rsid %in% top_markers$SNP)

table <- plot_ly(
  type = 'table',
  header = list(
    columnwidth = c(80,600),
    values = c(names(annotation_top)),
    align = c('left', rep('center', ncol(annotation_top))),
    line = list(width = 1, color = 'black'),
    fill = list(color = 'rgb(235, 100, 230)'),
    font = list(family = "Arial", size = 14, color = "white")
  ),
  cells = list(
    values = rbind(
      rownames(annotation_top),
      t(as.matrix(unname(annotation_top)))
    ),
    align = c('left', rep('center', ncol(annotation_top))),
    line = list(color = "black", width = 1),
    fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222, 249, 0.65)')),
    font = list(family = "Arial", size = 12, color = c("black")),
    domain = list(x=c(0,0.5), y=c(0,1))
  ))

# Plotting 
plotly::subplot(manhattan, volcano, qqplot, nrows = 3, margin = 0.1, titleX = T, titleY = T) %>%
  layout(title = "Manhattan, Volcano and QQ plots")#, height=(700))

main_plot <- plotly::subplot(manhattan, volcano, qqplot, nrows = 3, margin = 0.1, titleX = T, titleY = T, heights = c(0.4,0.3,0.2)) %>%
  layout(title = "Manhattan, Volcano and QQ plots")#, height=(700))

htmlwidgets::saveWidget(as.widget(main_plot), "multiqc_report.html")
htmlwidgets::saveWidget(as.widget(table), "annotated_table.html")

# Table annotated to display
write.table(annotation,  
            file      = 'annotation.csv',
            append    = FALSE, 
            quote     = FALSE, 
            sep       = ",",
            row.names = F,
            col.names = T)

# Make circular plot
CMplot(gwasResults, plot.type="c", r=1.6, cir.legend=TRUE,
        outward=TRUE, cir.legend.col="black", cir.chr.h=.1 ,chr.den.col="orange", file="jpg",
        memo="", dpi=300, chr.labels=seq(1,22))

#number of hits 
sprintf("Number of variants tested %i", nrow(gwas_input))
sprintf("Number of variants at suggestive significance %i", gwas_input %>% filter(P < 1e-5) %>% with(nrow(.)))
sprintf("Number of variants at suggestive significance %i", gwas_input %>% filter(P < 1e-8) %>% with(nrow(.)))