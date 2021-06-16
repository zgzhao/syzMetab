#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: mpath.R
# Description: handlng KEGG metabolic pathway (mpath) objects
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-05-31 12:12:17

#' Make metabolic pathway (keggPATH object)
#'
#' This function reads from local folder (d.path) a KEGG xml file, which will be downloaded if not exists. Make sure internet is accessible.
#' Refer to \code{\link{keggPATH}} for class information.
#' TODO: make pathway from local KGML file regardless of d.path
#' @title make metabolic pathway
#' @param object character, KEGG pathway ortholog (ko) index string such as ko00010, 00010 or ath00010
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return keggPATH object
#' @author ZG Zhao
#' @export
setGeneric("make_mpath", function(object, d.path="KEGG") standardGeneric("make_mpath"))
setMethod("make_mpath", "character", function(object, d.path){
    ko <- tolower(object[1])
    new("keggPATH", object, d.path)
})


#' Convert KOG names of a keggPATH object to real gene names
#'
#' In keggPATH objects, genes involved in reactions are generally represented by KEGG entry names (k genes or KOG).
#' @title General to organism-specific keggPATH conversion
#' @param object a keggPATH object
#' @param org abbreviation of an organism such as "ath"
#' @param d.path file path for \code{\link{KEGG_get}}
#' @return keggPATH object
#' @author ZG Zhao
#' @export
setGeneric("setorg_path", function(object, org, d.path="KEGG") standardGeneric("setorg_path"))
setMethod("setorg_path", c("keggPATH", "character"), function(object, org, d.path){
    if(Organism(object) != "ko") stop("No a general KEGG pathway (eg. ko00010)!")
    if( is.empty(org) ) return(object)
    org <- tolower(org)
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    ## filter reactions
    rns <- Reactions(object)
    rns <- lapply(rns, FUN=function(rr){
        gg <- intersect(rr$gene, kogs)
        gg <- if(length(gg) < 1) "absent" else unlist(gmap[gg])
        names(gg) <- NULL
        rr$gene <- .setGeneNames(gg)
        rr
    })
    ## ss <- sapply(rns, FUN=function(x) ! is.empty(x$gene))
    ## rns <- rns[ss]
    class(rns) <- "ReactionList"
    object@reactions <- rns
    ## filter entries
    ents <- object@entries
    ents <- lapply(ents, FUN=function(rr){
        if(! rr$type %in% "ortholog") return(rr)
        oo <- intersect(rr$name, kogs)
        oo <- if(length(oo) < 1) "absent" else unlist(gmap[oo])
        names(oo) <- NULL
        ox <- grepl("^K[0-9]{5}$", oo)
        rr$name <- sort(oo[!ox])
        rr$type <- "gene"
        rr
    })
    ## ss <- sapply(ents, FUN=function(x) ! is.empty(x$name))
    ## ents <- ents[ss]
    class(ents) <- "EntryList"
    object@entries <- ents
    object@pathInfo <- lapply(object@pathInfo, FUN=function(x) {
        gsub("ko([0-9]+)", paste0(org, "\\1"), x)
    })
    object@pathInfo$org <- org
    object
})
