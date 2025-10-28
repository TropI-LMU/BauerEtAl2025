## ---------------------------
##
## Script name: Get_Annotation.R
##
## Purpose of script: minor functions to retieve annotation from org.Mm/Hs.eg.Dbi
##
## Author: Olga Baranov
##
## Date Created: 2020-05-15
##
## Email: olga.baranov@posteo.de
##
## ---------------------------
##
## Notes: required libraries:
##                  org.Mm.eg.dbi
##                  org.Hs.eg.dbi
##                  dplyr
##
## ---------------------------


#' Retrieves info (i.e. column) by querying org.Xx.eg.dbi
#' @param genes list of gene IDs to retrieve info for; any ID present in org.Xx.eg.dbi are allowed
#' @param keytype type of ID used in the gene list; must tbe a column in org.Xx.eg.dbi
#' @param species either mouse or human
#' @param get.columns which columns
#' @param merge
#'

library(org.Hs.eg.db)
library(tidyverse)
library(biomaRt)


what.columns <- function(species){
  if(species == 'mouse' | species == 'mmu'){
    library(org.Mm.eg.db)
    keytypes(org.Mm.eg.db)
  } else if (species == 'human' | species == 'hsa') {
    library(org.Hs.eg.db)
    keytypes(org.Hs.eg.db)
  }
}

get.annotation <- function(genes, keytype, species, get.columns,  merge = 'none', merge.by = 'keytype'){
  del.ENS = FALSE
  if(keytype == 'SYMBOL'){
      library(limma)
      normed.genes = sapply(genes,function(x){alias2Symbol(x)[1]})
      genes = names(normed.genes) %>%
              sapply(function(x){if (is.na(normed.genes[[x]])) x else normed.genes[[x]][1]})
  }
  if(! 'ENSEMBL' %in% get.columns & keytype!='ENSEMBL'){get.columns <- c(get.columns,'ENSEMBL'); del.ENS <- TRUE}
  if(species == 'mouse' | species == 'mmu'){
    library(org.Mm.eg.db)
    ano.data <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, columns = get.columns, keytype = keytype)
  } else if (species == 'human' | species == 'hsa') {
    library(org.Hs.eg.db)
    ano.data <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, columns = get.columns, keytype = keytype)
  }
  # implement a warning for one to many mapping
  if(keytype == 'SYMBOL'){
      genes.dct = names(normed.genes)
      names(genes.dct) = genes
      ano.data$origSYMBOL = ano.data$SYMBOL %>%
                      sapply(function(x){genes.dct[x]}) %>% unname
  }
  if(merge.by == 'keytype'){merge.by = keytype}
  ano.data <- merge.annotation(ano.data, merge, merge.by)
  if(del.ENS){ano.data <- ano.data %>% dplyr::select(-'ENSEMBL')}
  return(ano.data)
}

get.description = function(genes, species){
  if (species == "mouse" | species == "mmu") {
    dataset = "mmusculus_gene_ensembl"
  } else if (species == "human" | species == "hsa") {
    dataset = "hsapiens_gene_ensembl"
  }
  ensembl = useEnsembl(biomart="ensembl", dataset=dataset)
  genedesc <- getBM(attributes = c("ensembl_gene_id", "description"), filters = "ensembl_gene_id",
                   values = genes, mart = ensembl)
  return(genedesc)
}


#' Retrieves info (i.e. column) by querying org.Xx.eg.dbi
#' @param genes list of gene IDs to retrieve info for; any ID present in org.Xx.eg.dbi are allowed
#' @param keytype type of ID used in the gene list; must tbe a column in org.Xx.eg.dbi
#' @param species either mouse or human
#' @param merge how to merge data in case of one to many mapping in database; available: none, first, all
#'
get.annotation.CHR <- function(genes, keytype, species, merge = 'none'){
  # outcommented regions can probably be deleted
  # if(keytype == 'ENTREZID'){
    # ano.data <- tibble(ENTREZID = genes)
  # } else
  if( species == 'mouse'){
    library(org.Mm.eg.db)
    ano.data <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, columns = c("ENTREZID"), keytype = keytype)
    get.chr <- function(x){sapply(x, function(y){ if(is.na(y)) NA else unlist(as.list(org.Mm.egCHR[y]))[1]})}
  } else if (keytype != 'ENTREZID' & species == 'human') {
    library(org.Hs.eg.db)
    ano.data <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, columns = c("ENTREZID"), keytype = keytype)
    get.chr <- function(x){sapply(x, function(y){ if(is.na(y)) NA else unlist(as.list(org.Hs.egCHR[y]))[1]})}
  }
  # ano.data <- merge.annotation(ano.data, merge, keytype)
  #helper fct to get the chromosome pos 
  ano.data <- tibble(ano.data) %>% mutate(chr = get.chr(ENTREZID))
  return(ano.data)
}

#' Retrieves info (i.e. column) by querying org.Xx.eg.dbi
#' @param genes list of gene IDs to retrieve info for; any ID present in org.Xx.eg.dbi are allowed
#' @param keytype type of ID used in the gene list; must tbe a column in org.Xx.eg.dbi
#' @param species either mouse or human
#' @param merge how to merge data in case of one to many mapping in database; available: none, first, all
#'
get.annotation.CHRPOS <- function(genes, keytype, species, merge = 'none'){
  if(keytype == 'ENTREZID'){
    ano.data <- tibble(ENTREZID = genes)
  } else if( species == 'mouse'){
    library(org.Mm.eg.db)
    ano.data <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, columns = c("ENTREZID"), keytype = keytype)
  } else if (keytype != 'ENTREZID' & species == 'human') {
    library(org.Hs.eg.db)
    ano.data <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, columns = c("ENTREZID"), keytype = keytype)
  }
  # ano.data <- merge.annotation(ano.data, merge, keytype)
  #helper fct to get the chromosome pos
  get.chr <- function(x){sapply(x, function(y){ if(is.na(y)) NA else unlist(as.list(org.Mm.egCHR[y]))[1]})}
  get.chr.start <- function(x){sapply(x, function(y){ if(is.na(y)) NA else unlist(as.list(org.Mm.egCHRLOC[y]))[1]})}
  get.chr.end <- function(x){sapply(x, function(y){ if(is.na(y)) NA else unlist(as.list(org.Mm.egCHRLOCEND[y]))[1]})}
  ano.data <- ano.data %>% mutate(chr = get.chr(ENTREZID), chr.pos.start = get.chr.start(ENTREZID), chr.pos.end = get.chr.end(ENTREZID) )
  return(ano.data)
}


## helper function
merge.annotation <- function(ano.data, merge, keytype){
  dups <- duplicated(ano.data[, 1])
  if(merge == 'first'){
    if(sum(dups >0 )){warning(paste0('multiple annotations found for: ', paste(unlist(ano.data[dups,1]), collapse=' ') ))}
    ano.data <- ano.data[! dups, ]
    # ano.data <- ano.data %>% group_by(!! sym(keytype)) %>% filter(row_number()==1) for a dplyrry way
  } else if (merge == 'all'){
    warning(paste0('WARNING: using merge = all joins annotations by semicolon, you might want to review the data before using the table for further analysis;',
                   'multiple annotations found for: ', paste(unlist(ano.data[dups,1]), collapse=' ')))
    merge.by.semicolon <- function(clmn){paste(clmn, collapse=';')}
    ano.data <- ano.data %>% group_by(!! sym(keytype)) %>% mutate_at(vars(-group_cols()), merge.by.semicolon) %>% ungroup() %>%  unique()
  } else if(merge == 'drop') {
    dup.genes <- unlist(ano.data[which(dups),1])
    if(sum(dups >0 )){warning(paste0('multiple annotations found for: ', paste(unlist(ano.data[dups,1]), collapse=' '), '; dropping the genes', collapse=' ' ))}
    ano.data <- ano.data  %>% filter(! ENSEMBL %in% dup.genes)
  } else {warning('~~~~~~~~~~~~~~~~~wrong filtering option~~~~~~~~~~~~~~~~~~~~~~~~~')}
  return(as_tibble(ano.data))
}


reduce.go = function(terms){
    library(GO.db)
    bpc = as.data.frame(GOBPCHILDREN)
    gobp = as_tibble(bpc[1:2]) %>%
        rename('go_id' = 'child', 'go_id.1' = 'parent')

    # filter which terms are in the parent column; those are leaf nodes
     parentinset = gobp %>% filter(parent %in% terms)

    # for the above check if a child is in the list, discard
    dropgo = parentinset %>%
        filter(child %in% terms) %>%
        pull(parent) %>%
        unique

    #discard those who have, keep the rest
    return(setdiff(terms, dropgo))
}




#' Converts EnsemblIDs to EntrezIDs using the above functions; can handle version notation
#' @param tab table with ensembl IDs as rownames
#' @param originalID the ID type in the rownames; (use what.columns to see options)
#' @param targetID the ID type to convert to; (use what.columns to see options)
#' 
replaceIDs = function(tab, originalID, targetID, keepOriginal = FALSE){
    if(originalID == 'ENSEMBL'){
          ensids = rownames(tab) %>%
            sapply(function(x){strsplit(x, '\\.')[[1]][1]})
    } else { ensids = rownames(tab) }
    tab = as_tibble(tab)
    tab['originalID'] = ensids

    ezid = tab %>%
            pull('originalID') %>%
            get.annotation( originalID,'human',targetID)
    ezid[is.na(ezid[targetID]),targetID] = ezid[is.na(ezid[targetID]),originalID]
    tmp = tab %>%
        inner_join(ezid, by = c('originalID'=originalID)) %>%
        drop_na(all_of(targetID)) %>%
        distinct(!!sym(targetID), .keep_all = TRUE) %>%
        column_to_rownames(targetID)
    tmp = (if (keepOriginal) tmp else (tmp %>% dplyr::select(-"originalID")))
    return(tmp)
}
