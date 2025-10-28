library(tidyverse)
library(vroom)
library(DESeq2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(pROC)
library(corrplot)

prefix <- ""
basedir <- "/path/to/Repo/"
outdir <- paste0(basedir, "/OutputData/")

source("TBSignatures.R")
source("Get_Annotation.R")


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


shapes <- c(
    # marker colors
    'Control' = 21,
    'Non-converter' = 23,
    'Reverter at EOT' = 22,
    'Recurrence after EOT' = 25
)




vsd <- readRDS(paste0(outdir, "/RData/vsdNormCounts_allreentered.Rdata"))
counts_norm <- assay(vsd)
rownames(counts_norm) = rownames(counts_norm) %>% sapply(function(x){str_split(x,'\\.')[[1]][1]}) %>% unname

colinfo_reenter <- vroom(paste0(basedir, "/FilteredData/colinfo_reentered.csv"),
    col_types = c(Sample_Name = "c")
) %>% dplyr::select(-'...1')

cd38 = read_xlsx(paste0(basedir,'/RawData/TAM_TB_CD38.xlsx'), sheet = 'BL_M12')

# extract the cd38 values for the respective time point 

group = colinfo_reenter %>% pull(type) %>% 
        sapply(function(x){ifelse(x == 'Reentered',1,0)}) %>% 
        unname
mamsscore <- MAMS_6_mean(counts_norm[, colinfo_reenter$Sample_Name], useensembl = TRUE)
roc_mams <- roc(group, mamsscore)
tss_mams = coords(roc_mams)
mams_cutof = tss_mams[which.max(tss_mams[,2:3] %>% rowSums()),1]
swescore <- sweeney3(counts_norm[, colinfo_reenter$Sample_Name], useensembl = TRUE)
roc_swe <- roc(group, swescore)
tss_swe = coords(roc_swe)
swe_cutof = tss_swe[which.max(tss_swe[,2:3] %>% rowSums()),1]
riskscore <- risk6(counts_norm[, colinfo_reenter$Sample_Name], useensembl = TRUE)
roc_risk <- roc(group, riskscore)
tss_risk = coords(roc_risk)
risk_cutof = tss_risk[which.max(tss_risk[,2:3] %>% rowSums()),1]


info_mini = colinfo_reenter[c('PersonID','R_group','Timepoint','Sample_Name')]


scores = merge(info_mini, DataFrame(mamsscore, swescore, riskscore) %>% as.data.frame() %>% rownames_to_column('Sample_Name'), by = 'Sample_Name')
tmpfun = function(x){
    str_replace(x,'M00','BL') %>%
    str_replace('M02','M2') %>%
    str_replace('M06','M6') 
}

scores['Month'] = scores$Timepoint %>% sapply(tmpfun) %>% unname
cd38 = cd38 %>% rename(PersID='PersonID')
fulltab = merge(cd38,scores, by = c('PersonID','Month'), how = '') %>% rename(RatioCD38posOfIFNypos = 'TAM_TB', mamsscore = 'MAMS6', swescore = 'Sweeney3', riskscore = 'RISK6') 


missing = anti_join(cd38,scores, by = c('PersonID','Month'))
# 16 samples that might be in (m 0,2,6) are missing


cm = cor(fulltab[c('TAM_TB','MAMS6','Sweeney3','RISK6')] %>% as.matrix)
ctab = cm %>% as.data.frame %>% rownames_to_column('varA') %>% 
    pivot_longer(
        -varA,
        names_to = "varB", 
        values_to = "correlation"
  ) 

ctab['textcolor'] = ctab$correlation > 0

cp = ggplot(ctab, aes(x= varA, y = varB)) + 
    geom_tile(aes(fill = correlation)) +
    geom_text(aes(label = round(correlation, 2), color = textcolor)) +
    scale_color_manual(values = c('white','white')) + 
    scale_fill_gradient2() +
    xlab('') + ylab('') +  guides(color = guide_none()) +
    theme_minimal()+ theme(legend.position = "bottom")


namemap = c('control' = "Control", 'non converter' = "Non-converter", 
            'failure at M6' = "Reverter at EOT", 'relapse' = "Recurrence after EOT")
fulltab['R_group'] = fulltab['R_group'] %>% sapply(function(x){namemap[x]})
fulltab['R_group'] = factor(fulltab['R_group'] %>% pull, 
    levels = c("Control", "Non-converter", "Reverter at EOT", "Recurrence after EOT"))


fulltab[c('PersonID','Timepoint','R_group','TAM_TB','MAMS6','Sweeney3','RISK6')] %>%
    as.data.frame %>%
    write_delim(paste0(outdir,'ScoreTab.csv'), delim = '\t')

## for the final plot:
# 3 line plots showing the correlation btw. TAM_TB and a signature, with the dots 
# shape coded for ctrl / relapse at 6 / relapse post 6
# panel showing the cor matrix as a plot
# cut off for TAM-TB 31.6 % (i.e. 0.316 on the x axis)
tamcut = 0.316



mp = ggscatter(fulltab, x= 'TAM_TB', y = 'MAMS6', 
        add = "reg.line", add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, shape = 'R_group', ,  fill = 'R_group') +  #shape = 'R_group',  fill = 'R_group'
        stat_cor(method = "pearson") +
        scale_shape_manual(values=shapes, name = '') + 
        scale_fill_manual(values = cols, name = '') + 
        xlab('TAM-TB') 

rp = ggscatter(fulltab, x= 'TAM_TB', y = 'RISK6', 
        add = "reg.line", add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, shape = 'R_group', ,  fill = 'R_group') +   #shape = 'R_group',  fill = 'R_group'
        stat_cor(method = "pearson") +
        scale_shape_manual(values=shapes, name = '') + 
        scale_fill_manual(values = cols, name = '') + 
        xlab('TAM-TB') 

sp = ggscatter(fulltab, x= 'TAM_TB', y = 'Sweeney3', 
        add = "reg.line", add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE,shape = 'R_group', ,  fill = 'R_group') +  
        # annotate("rect", fill = "grey", alpha = 0.3, 
        #     xmin = tamcut, xmax = Inf,
        #     ymin = -Inf, ymax = Inf) +
        # annotate("rect", fill = "grey", alpha = 0.3, 
        #     xmin = -Inf, xmax = Inf,
        #     ymin = -Inf, ymax = swe_cutof) +
        stat_cor(method = "pearson") +
        scale_shape_manual(values=shapes, name = '') + 
        scale_fill_manual(values = cols, name = '') + 
        xlab('TAM-TB') 
# does not make sense since we use all timepoints, i.e. controls can have values below threshold for early timepoints


legend_1 <- get_legend(cp)
legend_2 <- get_legend(sp)

fullplot = ggarrange(cp,sp,mp,rp, ncol = 2, nrow = 2, legend.grob = rbind(legend_1, legend_2), labels = c('A','B','C','D')) 
ggsave(paste0(outdir,'/Plots/correlationFigure.pdf'), plot = fullplot)




#################### combining Sweeney and TAM

fulltab['SweeneyTAM'] = fulltab['TAM_TB'] - fulltab['Sweeney3']/5

fulltab %>% head()
fulltab$R_time_mcoded = fulltab$R_time %>% sapply(function(x){str_replace(x, 'Mth ','M') %>% str_replace('0','')}) %>% unname()

subtab = fulltab[fulltab$Month == fulltab$R_time_mcoded,]
group = subtab %>% pull(R_group) %>% 
        sapply(function(x){ifelse(x == 'control',0,1)}) %>% 
        unname
roc_sam <- roc(group, subtab$SweeneyTAM)
roc_tam <- roc(group, subtab$TAM_TB)
roc_swe <- roc(group, subtab$Sweeney3)
plot(roc_swe)
dev.new()
roc.list = list(TAM_TB = roc_tam,Sweeney3 = roc_swe, combined = roc_sam)
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
    values = c("TAM_TB" = "#dd001c", "Sweeney3" = "#066af4", 'combined' = '#660066')) +
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

