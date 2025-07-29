#!/bin/env Rscript

# generate_cna_baf_plots.R
# This script generates copy number and B-allele frequency plots from a BAF bedgraph file and a CN ratio TSV file.
# It creates two PNG files: cna_plot.png and baf_plot.png.

script_version <- "1.0.0"

# suppress warnings and messages
suppressPackageStartupMessages({
  require(dplyr)
  require(ggplot2)
  require(cowplot)
  require(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)

if (args[1]=="--version"){
  cat(paste0(script_version,"\n"))
  quit(status = 0)
}

if (length(args) < 3) {
  stop("Usage: generate_cna_baf_plots.R <case_id> <baf_bedgraph> <cn_ratio_tsv>")
}

caseid <- args[1]
baf_bedgraph <- args[2]
cn_ratio_tsv <- args[3]

bafs <- read.table(baf_bedgraph,sep="\t",col.names = c("chrom","start","end","value")) %>%
  mutate(type=NA,category="baf",csyntax=NA,known_genes=NA,filters=NA) %>%
  select(chrom,start,category,type,value,csyntax,known_genes,filters)
bafs$filters <- NA
bafs$baf <- bafs$value
bafs$cnratio <- NA

# Read CN ratio TSV file
cn <- read.table(cn_ratio_tsv,skip=3,sep="\t",
                   col.names = c("chrom","start","end","name","value","improper_pairs"))

# make GRanges object for CN ratio
cn_gr <- GRanges(cn$chrom,IRanges(cn$start,cn$end))                   
cn_gr$score <- cn$value

# format cn
cn <- cn %>% mutate(start=round(start + (end-start)/2,0),type=NA,category="cn",csyntax=NA,known_genes=NA,filters=NA) %>%
  select(chrom,start,category,type,value,csyntax,known_genes,filters)
cn$filters <- NA
cn$baf <- NA
cn$cnratio <- cn$value

# For now, we will not perform the centering of the cnratio, as it depends on the SV calls from the JSON report.
# This can be added back later if needed.

dat <- rbind(cn,bafs) %>% 
             mutate(chrom=factor(chrom,levels=paste0("chr",c(1:22,"X","Y")),ordered=T),start=as.numeric(start)) %>%
             arrange(chrom,start) %>% mutate(value=ifelse(category %in% c("baf","cn"),NA,value))
  
dat$index <- 1:nrow(dat)

label_pos <- round(rowMeans(cbind(c(0,dat$index[which(dat$chrom != lag(dat$chrom))]),c(dat$index[which(dat$chrom != lag(dat$chrom))],max(dat$index)))),0)
labelDf <- data.frame(pos=label_pos,
                      label=as.character(dat[round(rowMeans(cbind(c(0,dat$index[which(dat$chrom != lag(dat$chrom))]),c(dat$index[which(dat$chrom != lag(dat$chrom))],max(dat$index)))),0),"chrom"]),
                      csyntax=NA,filters=NA)

p_cn <- ggplot(dat,aes(x=index,y=cnratio)) +
        geom_bin2d(bins=1000,na.rm=TRUE) +
        geom_vline(xintercept = dat$index[which(dat$chrom != lag(dat$chrom))],linetype=2,col="gray") +
        scale_x_continuous(limits=c(1,max(dat$index)),breaks=dat$index[which(dat$chrom != lag(dat$chrom))],labels = NULL,expand = expansion(mult=0.02),
                           sec.axis = dup_axis(breaks=labelDf$pos,labels = labelDf$label)) +
        scale_y_continuous(name="copy number ratio",
                           limits=c(min(-2,round(quantile(dat$cnratio,0.05,na.rm = T),0)),
                                    max(2,round(quantile(dat$cnratio,0.05,na.rm = T,0.95),0))),expand = expansion(mult=c(0.01,0.02))) +
        scale_fill_stepsn(values=c(0,0.0001,1),colors = c(NA,"darkblue","darkblue")) +
        ggtitle(paste0("Copy number for ",caseid)) +
        theme_classic() + theme(plot.title = element_text(hjust=0.5),legend.position = "none",axis.ticks.x = element_blank(),axis.text.y = element_text(color = "black",size=10),axis.text.x.bottom = element_blank(),
                                axis.title.x = element_blank(),axis.text.x = element_text(angle=90,color="black"),plot.margin = margin(t=20, l=10)) + coord_cartesian(clip = "off")

p_baf <- ggplot(dat,aes(x=index,y=baf)) +
        geom_bin2d(binwidth=c(2000,0.005),na.rm=TRUE) +
        geom_vline(xintercept = dat$index[which(dat$chrom != lag(dat$chrom))],linetype=2,col="gray") +
        scale_x_continuous(limits=c(1,max(dat$index)),breaks=dat$index[which(dat$chrom != lag(dat$chrom))],labels = NULL,expand = expansion(mult=0.02),
                           sec.axis=dup_axis()) +
        scale_y_continuous(name="B-allele frequency",
                           limits=c(0,1),expand = expansion(mult=0.01)) +
        scale_fill_stepsn(values=c(0,0.0001,1),colors = c(NA,"darkblue","darkblue")) +
        ggtitle(paste0("B-allele frequency for ",caseid)) +
        theme_classic() + theme(legend.position = "none",axis.ticks.x = element_blank(),axis.text.y = element_text(color = "black",size=10),plot.margin = margin(l=10),
                                axis.text.x = element_blank(),axis.title.x = element_blank()) + coord_cartesian(clip = "off")

ggsave("cna_plot.png", plot = p_cn, width = 11, height = 4, units = "in", dpi = 300)
ggsave("baf_plot.png", plot = p_baf, width = 11, height = 4, units = "in", dpi = 300) 