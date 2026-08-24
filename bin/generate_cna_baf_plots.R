#!/usr/bin/env Rscript

# generate_cna_baf_plots.R
# This script generates copy number and B-allele frequency plots from a BAF bedgraph file and a CN ratio TSV file.
# If inputs are missing, it generates placeholder PNGs to prevent pipeline failures.

script_version <- "1.0.0"

# Load optparse first to handle arguments before loading heavy libraries
suppressPackageStartupMessages(require(optparse))

# --- Argument Parsing Logic ---

option_list <- list(
  make_option(c("-i", "--id"), type="character", default=NULL, 
              help="Case ID (used for titles and filenames)", metavar="character"),
  make_option(c("-b", "--baf"), type="character", default=NULL, 
              help="Path to BAF bedgraph input file", metavar="file"),
  make_option(c("-c", "--cn"), type="character", default=NULL, 
              help="Path to CN ratio TSV input file", metavar="file"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE,
              help="Print the script version and exit")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# --- Check for Version Flag ---
if (opt$version) {
  cat(paste0(script_version, "\n"))
  quit(save="no", status=0)
}

# Helper function to check validity
inputs_valid <- function(opt) {
  if (is.null(opt$id) || is.null(opt$baf) || is.null(opt$cn)) return(FALSE)
  if (!file.exists(opt$baf)) return(FALSE)
  if (!file.exists(opt$cn)) return(FALSE)
  return(TRUE)
}

# --- Validation & Placeholder Logic ---

if (!inputs_valid(opt)) {
  # Determine a safe ID for filenames/logging (fallback to "placeholder" if ID is missing)
  safe_id <- if (!is.null(opt$id)) opt$id else "placeholder"
  
  # Print warning to stderr
  cat(sprintf("[GENERATE_CNA_BAF_PLOTS] Missing BAF/CNV inputs for %s; creating placeholder plots\n", safe_id), file=stderr())
  
  # Generate Placeholder CNA Plot
  png(paste0(safe_id, ".cna_plot.png"), width=1800, height=900, res=150)
  par(mar=c(0,0,0,0))
  plot.new(); text(0.5, 0.5, "CNA plot unavailable", cex=2)
  dev.off()
  
  # Generate Placeholder BAF Plot
  png(paste0(safe_id, ".baf_plot.png"), width=1800, height=900, res=150)
  par(mar=c(0,0,0,0))
  plot.new(); text(0.5, 0.5, "BAF plot unavailable", cex=2)
  dev.off()
  
  # Exit successfully (0) so pipeline does not crash
  quit(save="no", status=0)
}

# --- Normal Execution (Only loads if inputs are valid) ---

suppressPackageStartupMessages({
  require(dplyr)
  require(ggplot2)
  require(cowplot)
  require(GenomicRanges)
})

# Assign to variables
caseid <- opt$id
baf_bedgraph <- opt$baf
cn_ratio_tsv <- opt$cn

# Read BAF bedgraph file
bafs <- read.table(baf_bedgraph, sep="\t", col.names = c("chrom","start","end","value")) %>%
  mutate(type=NA, category="baf", csyntax=NA, known_genes=NA, filters=NA) %>%
  select(chrom, start, category, type, value, csyntax, known_genes, filters)
bafs$filters <- NA
bafs$baf <- bafs$value
bafs$cnratio <- NA

# Read CN ratio TSV file
cn <- read.table(cn_ratio_tsv, skip=3, sep="\t",
                   col.names = c("chrom","start","end","name","value","improper_pairs"))

# make GRanges object for CN ratio
cn_gr <- GRanges(cn$chrom, IRanges(cn$start, cn$end))                   
cn_gr$score <- cn$value

# format cn
cn <- cn %>% mutate(start=round(start + (end-start)/2, 0), type=NA, category="cn", csyntax=NA, known_genes=NA, filters=NA) %>%
  select(chrom, start, category, type, value, csyntax, known_genes, filters)
cn$filters <- NA
cn$baf <- NA
cn$cnratio <- cn$value

dat <- rbind(cn, bafs) %>% 
             mutate(chrom=factor(chrom, levels=paste0("chr", c(1:22, "X", "Y")), ordered=T), start=as.numeric(start)) %>%
             arrange(chrom, start) %>% mutate(value=ifelse(category %in% c("baf", "cn"), NA, value))
  
dat$index <- 1:nrow(dat)

label_pos <- round(rowMeans(cbind(c(0, dat$index[which(dat$chrom != lag(dat$chrom))]), c(dat$index[which(dat$chrom != lag(dat$chrom))], max(dat$index)))), 0)
labelDf <- data.frame(pos=label_pos,
                      label=as.character(dat[round(rowMeans(cbind(c(0, dat$index[which(dat$chrom != lag(dat$chrom))]), c(dat$index[which(dat$chrom != lag(dat$chrom))], max(dat$index)))), 0), "chrom"]),
                      csyntax=NA, filters=NA)

p_cn <- ggplot(dat, aes(x=index, y=cnratio)) +
        geom_bin2d(bins=1000, na.rm=TRUE) +
        geom_vline(xintercept = dat$index[which(dat$chrom != lag(dat$chrom))], linetype=2, col="gray") +
        scale_x_continuous(limits=c(1, max(dat$index)), breaks=dat$index[which(dat$chrom != lag(dat$chrom))], labels = NULL, expand = expansion(mult=0.02),
                           sec.axis = dup_axis(breaks=labelDf$pos, labels = labelDf$label)) +
         scale_y_continuous(name="copy number ratio",
                            limits=c(min(-2, round(quantile(dat$cnratio, 0.05, na.rm = TRUE), 0)),
                                     max(2, round(quantile(dat$cnratio, 0.95, na.rm = TRUE), 0))), expand = expansion(mult=c(0.01, 0.02))) +
        scale_fill_stepsn(values=c(0, 0.0001, 1), colors = c(NA, "darkblue", "darkblue")) +
        ggtitle(paste0("Copy number for ", caseid)) +
        theme_classic() + theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.ticks.x = element_blank(), axis.text.y = element_text(color = "black", size=10), axis.text.x.bottom = element_blank(),
                                axis.title.x = element_blank(), axis.text.x = element_text(angle=90, color="black"), plot.margin = margin(t=20, l=10)) + coord_cartesian(clip = "off")

p_baf <- ggplot(dat, aes(x=index, y=baf)) +
        geom_bin2d(binwidth=c(2000, 0.005), na.rm=TRUE) +
        geom_vline(xintercept = dat$index[which(dat$chrom != lag(dat$chrom))], linetype=2, col="gray") +
        scale_x_continuous(limits=c(1, max(dat$index)), breaks=dat$index[which(dat$chrom != lag(dat$chrom))], labels = NULL, expand = expansion(mult=0.02),
                           sec.axis=dup_axis()) +
        scale_y_continuous(name="B-allele frequency",
                           limits=c(0, 1), expand = expansion(mult=0.01)) +
        scale_fill_stepsn(values=c(0, 0.0001, 1), colors = c(NA, "darkblue", "darkblue")) +
        ggtitle(paste0("B-allele frequency for ", caseid)) +
        theme_classic() + theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.y = element_text(color = "black", size=10), plot.margin = margin(l=10),
                                axis.text.x = element_blank(), axis.title.x = element_blank()) + coord_cartesian(clip = "off")

ggsave(paste0(caseid, ".cna_plot.png"), plot = p_cn, width = 11, height = 4, units = "in", dpi = 300)
ggsave(paste0(caseid, ".baf_plot.png"), plot = p_baf, width = 11, height = 4, units = "in", dpi = 300)