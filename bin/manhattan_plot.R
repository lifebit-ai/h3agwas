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

# Test to update markers to UKB

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
                        "\nFunction: ",gwas_plot$Function,"\nMAF: ", gwas_plot$F_A, sprintf("\nP-value: %.2e", gwas_plot$P),
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

manhattan <- ggplotly(p, tooltip="text", width = 1400, height = 1200)

volcano <- ggplotly(p1, tooltip="text", width = 1400, height = 1200)

# qqplot <- qq(gwas_plot$P) # NOT INTERACTIVE

qqplot <- qqly(gwas_plot %>% filter(P > 0), main= "QQ plot", snp="SNP", gene = "Gene", width = 1400, height = 1200)

# Plotting 
sub1 <- plotly::subplot(volcano, qqplot)

#plotly::subplot(sub1,manhattan, table, nrows = 2, titleX = T, titleY = T)

main_plot <- plotly::subplot(sub1, manhattan, nrows = 2, titleX = T, titleY = T) %>%
  plotly::layout(title = "Manhattan, Volcano and QQ plots")#, height=(700))

# Annotation table
annotation_top <- annotation %>% left_join(., gwas_plot %>% select(SNP,P), by =c("rsid"="SNP")) %>%
  # select( Location,rsid,Reference,Alternative, P, Gene,Function, AA.Change, SIFT.Effect,Polyphen2_HDIV.Effect,MutationTaster.Effect, CADD.Score,gnomAD.AF.All, gnomAD.AF.NFE) %>%
  mutate_all(list(~replace_na(.,"")))  %>% dplyr::arrange(P) %>% distinct(rsid,.keep_all = T) %>%
  select(
    rsid,Chromosome,Location,Reference,Alternative,Function,Gene,Gene.function,GWAS.catalogue.hits,Pathway,P,AA.Change,Clinvar.Publication, Clinvar.Pathogenicity, Clinvar,
    SIFT.Effect,SIFT.Score,Polyphen2.HDIV.Effect,Polyphen2.HDIV.Score,Polyphen2.HVAR.Effect,Polyphen2.HVAR.Score,
    MutationTaster.Effect,MutationTaster.Score,MutationAssessor.Effect,MutationAssessor.Score,CADD.Score,Vertebrate.Conserved, Mammalian.Conversed,
    gnomAD.AF.All, gnomAD.AF.AFR, gnomAD.AF.AMR,gnomAD.AF.ASJ,gnomAD.AF.EAS,gnomAD.AF.FIN,gnomAD.AF.NFE,gnomAD.AF.OTH                
    
  )

colnames(annotation_top) <- c(
  "rsid","Chromosome","Location","Reference","Alternative","Function","Gene","Gene function","GWAS catalogue hits","Pathway","P","AA Change","Clinvar Publication", "Clinvar Pathogenicity", "Clinvar",
  "SIFT Effect","SIFT Score","Polyphen2 HDIV Effect","Polyphen2 HDIV Score","Polyphen2 HVAR Effect","Polyphen2 HVAR Score",
  "MutationTaster Effect","MutationTaster Score","MutationAssessor Effect","MutationAssessor Score","CADD Score","Vertebrate Conserved", "Mammalian Conversed",
  "gnomAD AF All", "gnomAD AF AFR", "gnomAD AF AMR","gnomAD AF ASJ","gnomAD AF EAS","gnomAD AF FIN","gnomAD AF NFE","gnomAD AF OTH"
)

write.csv(annotation_top,"annotation_top_markers.csv", quote = F, row.names = F)

htmlwidgets::saveWidget(as.widget(main_plot), "multiqc_report.html")

# Make circular plot
CMplot(gwasResults, plot.type="c", r=1.6, cir.legend=TRUE,
        outward=TRUE, cir.legend.col="black", cir.chr.h=.1 ,chr.den.col="orange", file="jpg",
        memo="", dpi=300, chr.labels=seq(1,22))

#number of hits 
sprintf("Number of variants tested %i", nrow(gwas_input))
sprintf("Number of variants at suggestive significance %i", gwas_input %>% filter(P < 1e-5) %>% with(nrow(.)))
sprintf("Number of variants at suggestive significance %i", gwas_input %>% filter(P < 1e-8) %>% with(nrow(.)))