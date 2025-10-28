library(vroom)
library(tidyverse)
library(DESeq2)
library(readxl)
library(vsn)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)

prefix <- ""
basedir <- "/path/po/Repo/"
outdir <- paste0(basedir, "/OutputData/")


source("TBSignatures.R")
source("Get_Annotation.R")

## functions

drawLineplot <- function(score_title, scorevector, colinfo) {
    tmp <- colinfo %>%
        filter(Sample_Name %in% names(scorevector)) %>%
        dplyr::select(c("Sample_Name", "PersonID", "Timepoint", "type", "R_time","R_group")) %>%
        arrange(by = !!!syms(c("PersonID", "R_time", "type", "Timepoint")))

    tmp["Score"] <- scorevector[tmp %>% pull("Sample_Name")]
    tmp["hybridType"] <- as.character(tmp["R_group"] %>% pull())
    tmp[tmp$hybridType == "non converter",'hybridType'] = "failure at M6"
    tmp["hybridType"] <- factor(tmp %>% pull("hybridType"), levels = c("control", "failure at M6", "relapse"))
    # insert plot HERE
    mean0 <- tmp %>%
        filter(Timepoint == 0) %>%
        pull(Score) %>%
        mean()
    mean6 <- tmp %>%
        filter(Timepoint == 6) %>%
        pull(Score) %>%
        mean()
    mean12 <- tmp %>%
        filter(Timepoint == 12) %>%
        pull(Score) %>%
        mean()
    tmp["Timepoint"] <- as.numeric(as.character(tmp %>% pull("Timepoint")))
    tmp["culture"] <- "culture converter"
    tmp[tmp$R_group == "non converter", "culture"] <- "non-converter"
    tmp[tmp["type"] == "Control", 'R_time'] = 6
    tmp["R_time"] = factor(tmp$R_time, levels = tmp["R_time"] %>% unique() %>% pull() %>% sort())
    # tmp$culture = factor(tmp$culture, levels = c('non-converter','converter'))
    hmp <- ggplot(tmp, aes(x = Timepoint, y = Score, group = PersonID)) +
        # geom_ribbon(aes(ymin=mean6, ymax=mean6), fill = "grey70") +
        geom_line(aes(color = culture)) +
        geom_point(aes(color = culture, shape = R_time)) +
        geom_hline(aes(yintercept = mean6, colour = "global average at month 6")) +
        # geom_hline(yintercept = mean0, color = 'red', linetype='dashed' ) +
        # geom_hline(yintercept = mean12, color = 'blue', linetype='dashed' ) +
        facet_grid(hybridType ~ .) +
        scale_x_continuous(breaks = c(0, 2, 6, 9, 12)) +
        ggtitle(score_title) +
        labs(x = 'timepoint', y = 'score', color = '', shape = 'relapse month') +
        theme_classic(base_size = 14) +
        scale_color_manual(breaks = c("culture converter", "non-converter", "global average at month 6"), values = c("black", "red", "grey")) +
        scale_shape_manual(values = c(16,17,15,2,3,8), 
                labels = c('cured and 6','9','12','18','24','36'))
    return(hmp)
}

sigScorePlot <- function(sig_name, sig_fun) {
    siglist <- sig_fun(counts_norm, useensembl = TRUE)
    plt <- drawLineplot(sig_name, siglist, colinfo_reenter)
    return(plt)
}

scorefcts <- c(
    risk6 = risk6,
    sweeney3 = sweeney3,
    MAMS_6 = MAMS_6
)

## data

vsd <- readRDS(paste0(outdir, "/RData/vsdNormCounts_allreentered.Rdata"))
counts_norm <- assay(vsd)
rownames(counts_norm) = rownames(counts_norm) %>% sapply(function(x){str_split(x,'\\.')[[1]][1]}) %>% unname

colinfo_reenter <- vroom(paste0(basedir, "/FilteredData/colinfo_reentered.csv"),
    col_types = c(Sample_Name = "c")
) %>% dplyr::select(-'...1')

colinfo_reenter$R_time <- colinfo_reenter$vismon

colinfo_reenter$Timepoint <- factor(colinfo_reenter$Timepoint, levels = unique(sort(colinfo_reenter$Timepoint)))
colinfo_reenter$type <- factor(sapply(colinfo_reenter$type, 
                                function(x){str_replace(x,' ', '')}), 
                            levels = c('Control','Reentered'))

colinfo_reenter$Timepoint <- colinfo_reenter$Timepoint %>%
    sapply(function(x) {
        as.numeric(str_split(x, "M")[[1]][2])
    }) %>%
    unname()


pltswe <- sigScorePlot('Sweeney3', scorefcts[["sweeney3"]])
pltrisk <- sigScorePlot('RISK6', scorefcts[["risk6"]])
pltmams <- sigScorePlot('MAMS6', scorefcts[["MAMS_6"]])

uberp <- ggarrange(pltswe + labs(tag = "A"), pltrisk + labs(tag = "B"),
                    pltmams + labs(tag = "C"), ncol = 3,
                    common.legend = TRUE, legend = "bottom")
ggsave(paste0(basedir, "/OutputData/Plots/FigureY.pdf"), plot = uberp,
                width = 15, height = 8)
# dev.off()
