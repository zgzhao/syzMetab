#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: a_classes.R
# Description:
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-02-12 12:03:13

#' @name class_virtual
#' @title virtual classes
#' @aliases igraph mgraph KEGGdata
#' @description Virtual classes designed in mgraph package.
#' @details There following classes are set as virutal:
#' - igraph
#' - mgraph
#' - KEGGdata: actually list
#' @author zhao
NULL

## virtual classes
setClassUnion("igraph", "list")
setClassUnion("mgraph", "list")
setClassUnion("KEGGdata", "list")

#' @name class_KEGGmeta
#' @title S4 class: KEGGmeta
#' @aliases KEGGmeta
#' @seealso \code{\link{KEGGdata}}, \code{\link{MReactions}}
#' @description "KEGGmeta" is S4 class designed for holding KEGG pathway information: entries, reactions, relations, graphics and other general information (pathInfo).
#' @details Typically, a KEGGmeta object can be obtained by \code{\link{kmeta_from_ko}} function, which will download and parse the KEGG xml file.
#' @author zhao
NULL

setClass("KEGGmeta",
         slots=c(pathInfo="list",
                 entries="KEGGdata",
                 reactions="KEGGdata",
                 relations="KEGGdata",
                 graphics="KEGGdata"))

setMethod("initialize", "KEGGmeta", function(.Object, ko, d.path) {
    if(missing(ko) || is.empty(ko)) {
        return(.Object)
    }
    ko <- tolower(ko)
    f.path <- KEGG_get(ko, d.path, f.type="kgml")
    kdata <- read_xml(f.path)
    kinfo <- xml_children(kdata)
    if(length(kinfo) > 0) {
        entries <- .parseEntryList(kinfo)
        .Object@entries   <- entries
        .Object@relations <- .parseRelationList(kinfo)
        .Object@reactions <- .parseReactionList(kinfo, entries)
        .Object@graphics  <- .parseGraphicsList(kinfo)
    }
    .Object@pathInfo <- as.list(xml_attrs(kdata))
    return(.Object)
})


setMethod("show", "KEGGmeta", function(object) {
    cat("KEGGmeta object: S4 class holding KEGG pathway meta data\n")
    cat("\tName: ", object@pathInfo$name, "\n")
    cat("\tTitle: ", object@pathInfo$title, "\n")
    cat("\tOrganism: ", object@pathInfo$org,"\n")
    cat("\tReactions: ", length(object@reactions), "\n")
    cat("\tCompounds: ", length(getCPDs(object)), "\n")
    cat("\tGenes/Orthologs: ", length(getGenes(object)), "\n")
})

#' Refer to \code{\link{KEGGmeta}} for class information.
#'
#' This function reads from local folder (d.path) a KEGG xml file, which will be downloaded if not exists. Make sure internet is accessible.
#' @title Generate KEGGmeta object
#' @param ko character, KEGG pathway ortholog (ko) index string such as ko00010, 00010 or ath00010
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return KEGGmeta object
#' @author ZG Zhao
#' @export
kmeta_from_ko <- function(ko, d.path="KEGG"){
    ko <- tolower(ko[1])
    new("KEGGmeta", ko, d.path)
}


#' Convert KOG names of a KEGGmeta object to real gene names
#'
#' In KEGGmeta objects, genes involved in reactions are generally represented by KEGG entry names (k genes or KOG).
#' @title General to organism-specific KEGGmeta conversion
#' @param kinfo a KEGGmeta object
#' @param org abbreviation of an organism such as "ath"
#' @param d.path file path for \code{\link{KEGG_get}}
#' @return KEGGmeta object
#' @author ZG Zhao
#' @export
kmeta_x_org <- function(kinfo, org=NA, d.path="KEGG"){
    if(kinfo@pathInfo$org != "ko") stop("No a general KEGG pathway (eg. ko00010)!")
    if( is.empty(org) ) return(kinfo)
    org <- tolower(org)
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    ## filter reactions
    rns <- kinfo@reactions
    rns <- lapply(rns, FUN=function(rr){
        gg <- intersect(rr$gene, kogs)
        gg <- if(length(gg) < 1) NA else unlist(gmap[gg])
        names(gg) <- NULL
        rr$gene <- gg
        rr
    })
    ss <- sapply(rns, FUN=function(x) ! is.empty(x$gene))
    kinfo@reactions <- rns[ss]
    ## if(reaction.only) return(kinfo)

    ## filter entries
    ents <- kinfo@entries
    ents <- lapply(ents, FUN=function(rr){
        if(! rr$type %in% "ortholog") return(rr)
        gg <- intersect(rr$name, kogs)
        gg <- if(length(gg) < 1) NA else unlist(gmap[gg])
        names(gg) <- NULL
        rr$name <- gg
        rr$type <- "gene"
        rr
    })
    ss <- sapply(ents, FUN=function(x) ! is.empty(x$name))
    kinfo@entries <- ents[ss]
    kinfo@pathInfo <- lapply(kinfo@pathInfo, FUN=function(x) {
        gsub("ko([0-9]+)", paste0(org, "\\1"), x)
    })
    kinfo@pathInfo$org <- org
    kinfo
}


#' @name class_MReactions
#' @title S4 class: MReactions
#' @description MReactions: S4 class holding reaction infos.
#' @details Three slots: reaction, organism and alias.
#' - reaction: list of reaction info (substrate, product and genes)
#' - organism: KEGG organism identifier
#' - alias: named list of chemical aliases. Set when merging chemicals.
#' @author ZG Zhao
NULL

setClass("MReactions",
         slots=c(reaction="list",
                 organism="character",
                 alias="list"))
setMethod("initialize", "MReactions", function(.Object, reactions, organism="ko"){
    .Object@reaction <- reactions
    .Object@organism <- organism
    .Object@alias <- list()
    return(.Object)
})

setMethod("show", "MReactions", function(object) {
    cat("Reactions: S4 object holding KEGG reaction data\n")
    cat("\tOrganism: ", object@organism,"\n")
    cat("\tReactions: ", length(getReactions(object)), "\n")
    cat("\tCompounds: ", length(getCPDs(object)), "\n")
})

names.MReactions <- function(x){
    names(getReactions(x))
}
length.MReactions <- function(x){
    length(getReactions(x))
}
