library(vroom)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library("ggnewscale")
library(pROC)
library(gridExtra)


prefix <- ""
basedir <- "/path/to/Repo/"
outdir <- paste0(basedir, "/OutputData/")

source("TBSignatures.R")
source("Get_Annotation.R")


# data import and formating
# this table is created in the heatmapRelapseTimecourse.rmd script

vsd <- readRDS(paste0(outdir, "/RData/vsdNormCounts_allreentered.Rdata"))
counts_norm <- assay(vsd)
rownames(counts_norm) = rownames(counts_norm) %>% sapply(function(x){str_split(x,'\\.')[[1]][1]}) %>% unname

colinfo_reenter <- vroom(paste0(basedir, "/FilteredData/colinfo_reentered.csv"),
    col_types = c(Sample_Name = "c")
) %>% dplyr::select(-'...1')


# the plot itself
### violin plots

scorefcts <- c(
    risk6 = risk6,
    sweeney3 = sweeney3,
    MAMS_6_mean = MAMS_6_mean
)

# in new bernadettes tab R_time is missing

colinfo_reenter$Timepoint = colinfo_reenter$Timepoint %>% 
    sapply(function(x) {
        as.numeric(str_split(x, "M")[[1]][2])
    }) %>% unname

colinfo_reenter$R_time = colinfo_reenter$vismon
# some ctrls have very late timepoints of which we dont have the sequencing
# for them we can grab the latest available
for (ctrl in colinfo_reenter %>% dplyr::filter(Recurrence_Control == "C") %>% pull(PersonID) %>% unique()) {
    pats = colinfo_reenter$PersonID == ctrl
    subtab = colinfo_reenter %>% dplyr::filter(PersonID == ctrl)
    reltp = subtab$vismon[1]
    latetp = max(subtab$Timepoint)
    colinfo_reenter[pats, "R_time"] = latetp
}



colinfo_reenter['colgroup'] = ''
colinfo_reenter[colinfo_reenter$type == "Control", "colgroup"] = "control"
colinfo_reenter[colinfo_reenter$type == "Reentered" & colinfo_reenter$R_time == 6, 
                "colgroup"] = "failEOT"
colinfo_reenter[colinfo_reenter$type == "Reentered" & colinfo_reenter$R_time > 6, 
                "colgroup"] = "failafterEOT"
# this line has to go last
colinfo_reenter[colinfo_reenter$type == "Reentered" & colinfo_reenter$R_group == 'non converter', 
                "colgroup"] = "nonconverter"
colinfo_reenter["colgroup"] = factor(colinfo_reenter$colgroup, 
            levels = c('control','nonconverter','failEOT','failafterEOT'))
colinfo_reenter["type"] = factor(colinfo_reenter$type %>% 
            sapply(function(x){ifelse(x == 'Control','cured', 'recurrent')}) %>%
            unname, 
        levels = c('cured','recurrent'))


cols <- c(
    # marker colors
    'Control' = "#006633", # dark green
    "Non-converter" = "#ffee00", # yellow
    "Reverter at EOT" = "#ec3f0a", # light red
    "Recurrence after EOT" = "#7a1ae7", # purple
    # box colors
    'recurrent' = "#e42e46", # red
    'cured' = "#006633" # dark green
)


vioplot_altcolor <- function(counts, colinfo, scorefct, signame = "") {
    sub <- colinfo %>%
        dplyr::select(c(PersonID,R_time, type, Timepoint, Sample_Name, colgroup)) %>%
        filter(Timepoint <= R_time)
    siglist <- scorefct(counts[, sub$Sample_Name], useensembl = TRUE)
    tmp <- tibble(Sample_Name = names(siglist), score = siglist %>% unname())
    colnames(tmp) <- c("Sample_Name", signame)

    # sub$Timepoint <- factor(sub$Timepoint %>% unname(),
    #     levels = c('BL', 'M2', 'M6', 'M9', 'M12')
    # )
    tmp <- tmp %>% left_join(sub, by = "Sample_Name")
    tmp[tmp["type"] == "cured", "R_time"] = 6
    tmp["R_time"] <- factor(tmp$R_time, levels = c(6, 9, 12, 18, 24,36))
    namemap = c('converter' = "Control", 'nonconverter' = "Non-converter", 
                'failEOT' = "Reverter at EOT", 'failafterEOT' = "Recurrence after EOT")
    tmp['colgroup'] = tmp['colgroup'] %>% sapply(function(x){namemap[x]})
    tmp['colgroup'] = factor(tmp['colgroup'] %>% pull, 
        levels = c("Control", "Non-converter", "Reverter at EOT", "Recurrence after EOT"))
    tmp["Timepoint"] = tmp["Timepoint"] %>% 
        sapply(function(x){paste0('M',x) %>% str_replace('M0','BL')})

    tmp["Timepoint"] = factor(tmp["Timepoint"] %>% pull(), levels = c('BL','M2','M6','M9','M12'))
    
    p <- ggplot(tmp, aes_string(x = "Timepoint", y = signame, fill = "type")) +
        geom_boxplot(aes(fill = type), alpha = 0.3, outlier.shape = NA) +
        scale_fill_manual(
            name = "Group",
            values = cols[5:6],
            labels = c("C", "R")
        ) +
        guides(fill = guide_legend(nrow = 2)) +
        new_scale("fill") +
        geom_jitter(aes(group = type, shape = colgroup, fill = colgroup), # , color = colgroup
            position = position_jitterdodge(0.9), size = 1, alpha = 0.8
        ) +
        scale_shape_manual(
            values = c(21, 23, 22, 25),
            labels = c("Control", "Non-converter", "Reverter at EOT", "Recurrence after EOT")
        ) +
        labs(y = "score", x = "Time", color = "", shape = "Subgroup", fill = "Subgroup") +
        ggtitle(signame) +
        scale_fill_manual(
            values = cols[1:4],
            labels = c("Control", "Non-converter", "Reverter at EOT", "Recurrence after EOT")
        ) +
        guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
        facet_grid(cols = vars(colgroup))

    return(p + theme_classic() + theme(legend.position = "right"))
}


pmams <- vioplot_altcolor(counts_norm, colinfo_reenter, MAMS_6_mean, signame = "MAMS6")
# ggsave("/home/obaranov/projects/TBSequel/Relapses/OutputData/Plots/Violins/MAMS_6_mean.pdf", plot = pmams)
pswe <- vioplot_altcolor(counts_norm, colinfo_reenter, sweeney3, signame = "Sweeney3")
# ggsave("/home/obaranov/projects/TBSequel/Relapses/OutputData/Plots/Violins/sweeney3.pdf", plot = pswe)
prisk <- vioplot_altcolor(counts_norm, colinfo_reenter, risk6, signame = "RISK6")
# ggsave("/home/obaranov/projects/TBSequel/Relapses/OutputData/Plots/Violins/risk6.pdf", plot = prisk)





### ROCs ###################
############################
# prediction is the relapse / non relapse group 
# at timepoint of relapse
# m6
sub = colinfo_reenter %>% dplyr::filter(R_time == 6) %>% dplyr::filter(R_time == Timepoint)
sams = sub %>% pull(Sample_Name)
group = sub %>% pull(type) %>% 
        sapply(function(x){ifelse(x == 'recurrent',1,0)}) %>% 
        unname
mamsscore <- MAMS_6_mean(counts_norm[, sub$Sample_Name], useensembl = TRUE) %>% unname
roc_mams <- roc(group, mamsscore)
swescore <- sweeney3(counts_norm[, sub$Sample_Name], useensembl = TRUE)
roc_swe <- roc(group, swescore)
riskscore <- risk6(counts_norm[, sub$Sample_Name], useensembl = TRUE)
roc_risk <- roc(group, riskscore)


roc.list = list(MAMS6 = roc_mams,Sweeney3 = roc_swe, RISK6 = roc_risk)
ci.list <- lapply( roc.list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

p <- ggroc(roc.list) + theme_minimal() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal()

for(i in 1:3) {
  p <- p + geom_ribbon(
    data = dat.ci.list[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = i + 1,
    alpha = 0.2,
    inherit.aes = F) 
  } 

roc6 = p + scale_color_manual(name = "Signature", 
    values = c("MAMS6" = "#dd001c", "Sweeney3" = "#066af4", 'RISK6' = '#00b571')) +
    xlab("1 - specificity") + ylab("sensitivity") +
    coord_fixed() + theme_classic() +
    guides(color = guide_legend(ncol = 1))





# m6plus
sub = colinfo_reenter %>% dplyr::filter(R_time > 6) %>% dplyr::filter(R_time == Timepoint)
sams = sub %>% pull(Sample_Name)
group = sub %>% pull(type) %>% 
        sapply(function(x){ifelse(x == 'recurrent',1,0)}) %>% 
        unname
mamsscore <- MAMS_6_mean(counts_norm[, sub$Sample_Name], useensembl = TRUE) %>% unname
roc_mams <- roc(group, mamsscore)
swescore <- sweeney3(counts_norm[, sub$Sample_Name], useensembl = TRUE)
roc_swe <- roc(group, swescore)
riskscore <- risk6(counts_norm[, sub$Sample_Name], useensembl = TRUE)
roc_risk <- roc(group, riskscore)


roc.list = list(MAMS6 = roc_mams,Sweeney3 = roc_swe, RISK6 = roc_risk)
ci.list <- lapply( roc.list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

p <- ggroc(roc.list) + theme_minimal() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal()

for(i in 1:3) {
  p <- p + geom_ribbon(
    data = dat.ci.list[[i]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = i + 1,
    alpha = 0.2,
    inherit.aes = F) 
  } 

rocafter = p + scale_color_manual(name = "Signature", 
    values = c("MAMS6" = "#dd001c", "Sweeney3" = "#066af4", 'RISK6' = '#00b571')) +
    xlab("1 - specificity") + ylab("sensitivity") +
    coord_fixed() + theme_classic()  +
    guides(color = guide_legend(ncol = 1))


miep = grid.arrange(grobs = list(pmams + theme(legend.position="none"), 
                            pswe + theme(legend.position="none"), 
                            prisk+ theme(legend.position="none") ), 
                            nrow = 3)

legendscat = get_legend(pswe)
legendrocs = get_legend(roc6)

moep = grid.arrange(grobs = list(roc6 + theme(legend.position="none"), 
                            rocafter+ theme(legend.position="none"),
                            grid.arrange(grobs = list(legendscat, legendrocs), ncol = 2) ), 
                            nrow = 3)



troet = arrangeGrob(grobs = list(miep,moep), 
            layout_matrix = rbind(c(1,1,1,2,2),
                        c(1,1,1,2,2))
            ) 

ggsave(paste0(outdir, '/Plots/FigureX_new.pdf'), plot = troet, 
    width = 24,
    height = 18,
    units = "cm")
