## ---------------------------
##
## Script name: TBSignatures.R
##
## Purpose of script: functions to calculate parsimonous signatures as listed in Mendelson paper
##
## Author: Olga Baranov
##
## Date Created: 2023-05-25
##
## Email: olga.baranov@posteo.de
##
## ---------------------------
##
## Notes: required libraries:
##                  dplyr
##
## ---------------------------




#' All signatures were designed for raw CT values 
#' @params counts - a count table
#' @params useensembl - if set to TRUE, ensembl IDs are used, otherwise - symbols
#'


library(robustbase)

nthroot = function(x,n) {
  (abs(x)^(1/n))*sign(x)
}


francisco2 = function(counts, useensembl = FALSE){
    if(useensembl){
        counts['ENSG00000127528',] - counts['ENSG00000154451',]
    } else {
        counts['KLF2',] - counts['GBP5',] 
    }
}

maertzdorf4 = function(counts, useensembl = FALSE){
    if(useensembl){
        counts['ENSG00000117318',] - c( (counts['ENSG00000117228',] + counts['ENSG00000142089',] + counts['ENSG00000174944',]) / 3 )
    } else {
       counts['ID3',] - c( (counts['GBP1',] + counts['IFITM3',] + counts['P2RY14',]) / 3 )
    }
}

risk6 <- function(counts, useensembl = FALSE) {
    if (useensembl) {
        nthroot(counts["ENSG00000162645", ] * counts["ENSG00000198019", ] * counts["ENSG00000149131", ], 3) - nthroot(counts["ENSG00000128159", ] * counts["ENSG00000099899", ] * counts["ENSG00000100445", ], 3)
    } else {
        #  maybe FCGR1BP (with P, as it is a pseudogene)
        tryCatch(
            {nthroot(counts["GBP2", ] * counts["FCGR1BP", ] * counts["SERPING1", ], 3) - nthroot(counts["TUBGCP6", ] * counts["TRMT2A", ] * counts["SDR39U1", ], 3)},
        error = function(cond) { # switch to old gene name for FCGR1BP
            nthroot(counts["GBP2", ] * counts["FCGR1B", ] * counts["SERPING1", ], 3) - nthroot(counts["TUBGCP6", ] * counts["TRMT2A", ] * counts["SDR39U1", ], 3)}
    )}
}


roe1 = function(counts, useensembl = FALSE){
    if(useensembl){
        counts['ENSG00000168062',] -  nthroot(counts['ENSG00000139644',] * counts['ENSG00000070831',] * counts['ENSG00000105698',] * counts['ENSG00000115091',],4)
    } else {
        counts['BATF2',] - nthroot(counts['TMBIM6',] * counts['CDC42',] * counts['USF2',] * counts['ACTR3',],4)
    }
}

roe3 = function(counts, useensembl = FALSE){
    if(useensembl){
        ( (counts['ENSG00000074660',] + counts['ENSG00000154451',] + counts['ENSG00000168062',]) / 3
            - nthroot(counts['ENSG00000139644',] * counts['ENSG00000070831',] * counts['ENSG00000105698',] * counts['ENSG00000115091',],4) )
    } else {
        ((counts['SCARF1',] + counts['GBP5',] + counts['BATF2',]) / 3
            -  nthroot(counts['TMBIM6',] * counts['CDC42',] * counts['USF2',] * counts['ACTR3',],4) )
    }
}

sweeney3 = function(counts, useensembl = FALSE){
    if(useensembl){
        counts['ENSG00000127528',] - c( (counts['ENSG00000154451',] + counts['ENSG00000108861',]) / 2 )
    } else {
        counts['KLF2',] - c((counts['GBP5',] + counts['DUSP3',]) / 2 )
    }
}

MAMS_6_mean <- function(counts, useensembl = FALSE) {
    if (useensembl) {
        genes = c("ENSG00000163568", "ENSG00000198019", "ENSG00000120217", "ENSG00000117228", "ENSG00000082014", "ENSG00000119686")
    } else {
        genes = intersect(rownames(counts), c("AIM2", "FCGR1BP", "FCGR1B", "CD274", "GBP1", "SMARCD3", "FLVCR2"))
    }
    return(colMeans(counts[genes, ]))
}


MAMS_6p3_mean <- function(counts, useensembl = FALSE) {
    if (useensembl) {
        genes <- c("ENSG00000163568", "ENSG00000198019", "ENSG00000120217", "ENSG00000117228", "ENSG00000082014", 
                    "ENSG00000119686", "ENSG00000166710", "ENSG00000165704", "ENSG00000026025")
    } else {
        # "FCGR1B" probably has "FCGR1BP" as a more common name
        genes <- intersect(rownames(counts), c("AIM2", "FCGR1BP", "FCGR1B", "CD274", "GBP1", "SMARCD3", "FLVCR2", "B2M", "HPRT1", "VIM"))
    }
    return(colMeans(counts[genes, ]))
}


MAMS_6_med <- function(counts, useensembl = FALSE) {
    if (useensembl) {
        genes <- c("ENSG00000163568", "ENSG00000198019", "ENSG00000120217", "ENSG00000117228", "ENSG00000082014", "ENSG00000119686")
    } else {
        genes <- c("AIM2", "FCGR1BP", "CD274", "GBP1", "SMARCD3", "FLVCR2")
    }
    return(colMedians(counts[genes, ]))
}


MAMS_6p3_med <- function(counts, useensembl = FALSE) {
    if (useensembl) {
        genes <- c(
            "ENSG00000163568", "ENSG00000198019", "ENSG00000120217", "ENSG00000117228", "ENSG00000082014",
            "ENSG00000119686", "ENSG00000166710", "ENSG00000165704", "ENSG00000026025"
        )
    } else {
        # "FCGR1B" probably has "FCGR1BP" as a more common name
        genes <- c("AIM2", "FCGR1BP", "CD274", "GBP1", "SMARCD3", "FLVCR2", "B2M", "HPRT1", "VIM")
    }
    return(colMedians(counts[genes, ]))
}

# technically "FCGR1BP" and "FCGR1B" are the same gene, latter being the old naming convention
# but it is very unlikely that a data set will have both, so its a quick fix to just look for both

MAMS6p3genes = list(ensembl = c(
    "ENSG00000163568", "ENSG00000198019", "ENSG00000120217", "ENSG00000117228", "ENSG00000082014",
    "ENSG00000119686", "ENSG00000166710", "ENSG00000165704", "ENSG00000026025"
    ), symbol = c("AIM2", "FCGR1BP","FCGR1B", "CD274", "GBP1", "SMARCD3", "FLVCR2","B2M","HPRT1","VIM"))

MAMS6genes <- list(ensembl = c("ENSG00000163568", "ENSG00000198019", "ENSG00000120217", "ENSG00000117228",
     "ENSG00000082014", "ENSG00000119686"), 
     symbol = c("AIM2", "FCGR1BP", "FCGR1B","CD274", "GBP1", "SMARCD3", "FLVCR2"))

sweeney3genes = list(ensembl = c('ENSG00000127528','ENSG00000154451','ENSG00000108861'), 
    symbol = c('KLF2','GBP5','DUSP3'))

roe1genes = list( ensembl = c('ENSG00000168062','ENSG00000139644','ENSG00000070831','ENSG00000105698','ENSG00000115091'), symbol =c('BATF2','TMBIM6','CDC42','USF2','ACTR3') )

risk6genes = list(ensembl = c("ENSG00000162645", "ENSG00000198019", "ENSG00000149131", "ENSG00000128159", "ENSG00000099899", "ENSG00000100445"), symbol = c("GBP2", "FCGR1BP", "FCGR1B", "SERPING1", "TUBGCP6", "TRMT2A", "SDR39U1"))