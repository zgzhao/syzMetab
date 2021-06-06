#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: mgraph.R
# Description: handling metabolic graph (mgraph) objects
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-05-31 12:10:39

#' @name make_mgraph
#' @title make graph from various data
#' @aliases make_rgraph make_ggraph
#' @description Make metabolic graph from various type of data: keggPATH, ReactionSet, ReactionList, list or KOs (KEGG pathway identifiers).
#' @details Metabolic network can be represented by metabolite, reaction or gene graph.
#' @usage
#' - make_mgraph(object, ...)
#' - make_rgraph(object, ...)
#' - make_ggraph(object, ...)
#' @param object keggPATH, ReactionSet ReactionList, character vector or list.
#' @param org character, organism/species indentifier, parameter for "object" without organism info: ReactionList and list
#' @param d.path character, refer to \code{\link{KEGG_get}} for detail. Used when "object" is a KO vecter (KOs).
#' @author zhao
NULL

#' @export
setGeneric("make_mgraph", function(object, ...) standardGeneric("make_mgraph"))
setMethod("make_mgraph", "ReactionSet", function(object){
    ## basic method called by other methods!
    vv <- Compounds(object)
    g <- make_empty_graph(n=length(vv))
    V(g)$name <- vv
    rlist <- Reactions(object)
    for(ndx in names(rlist)) {
        rx <- rlist[[ndx]]
        genes <- rx[["gene"]]
        ss <- rx[["substrate"]]
        pp <- rx[["product"]]
        g <- add.edges(g, expand.grid(ss, pp), reaction=ndx)
    }
    attr(g, "reactions") <- object
    class(g) <- c("mgraph", class(g))
    g
})
setMethod("make_mgraph", "keggPATH", function(object){
    rsobj <- as_rset(object)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "character", function(object, d.path="KEGG"){
    rsobj <- rset_from_kos(object, d.path)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "ReactionList", function(object, org){
    rsobj <- as_rset(object, org)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "list", function(object, org){
    rvalid <- sapply(object, FUN=function(x){
        all(c("substrate", "product", "gene", "reversible") %in% names(x))
    })
    if(! all(rvalid)) stop("Invalid reaction list.")
    rsobj <- as_rset(rtns, org)
    make_mgraph(rsobj)
})

#' Normalize KEGG generic pathway to species specific pathway.
#'
#' Adjust and filter reactions, and rebuild graph object from reactions list.
#' @title general to organism-specific mgraph conversion
#' @param g mgraph object
#' @param org character
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return mgraph with organism set
#' @author ZG Zhao
#' @export
mgraph_x_org <- function(g, org, d.path = "KEGG"){
    if(! is.mgraph(g)) stop("Not a metabolic graph.")
    if(Organism(g) != "ko") stop("Not a generic metabolic graph!")

    rtns <- Reactions(g)
    org <- org[1]
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    rtns <- lapply(rtns, FUN=function(rr){
        gg <- intersect(rr$gene, kogs)
        gg <- if(length(gg) < 1) NA else unlist(gmap[gg])
        names(gg) <- NULL
        rr$gene <- gg
        rr$reversible <- FALSE
        rr
    })
    ss <- sapply(rtns, FUN=function(x) ! is.empty(x$gene))
    ## ensure returning result is a graph!
    if(sum(ss) < 1) return(.emptyGraph(g))

    rtns <- rtns[ss]
    robj <- as_rset(rtns, org)
    gx <- make_mgraph(robj)
    if(is.chemset(g)) {
        s <- intersect(Substrates(g), vnames(gx))
        p <- intersect(Products(g), vnames(gx))
        if(is.empty(s) || is.empty(g))
            gx <- .emptyGraph(gx)
        else
            gx <- xgraph_setchems(gx, s, p, clean=FALSE)
    }
    gx
}

#' Append additional reactions to KEGG pathway.
#'
#' TODO: may not work
#' Reactions without KEGG index cannot be extrated from KGML files. This function help you add these reactions manually.
#' @title append reactions to mgraph manually
#' @param g graphMET object
#' @param df a data frame with at least four columns which supply "from", "to", "reversible", "genes" infos for the reactions. File path to a csv file containing such infos is also acceptable.
#' @return
#' @author ZG Zhao
#' @export
mgraph_append <- function(g, df) {
    if(! is.data.frame(df) && file.exists(df)) df <- read.csv(df)
    if(ncol(df) < 4) stop("Data should have at least four columns: from, to, reversible and genes.")
    if(all(c("from", "to", "reversible", "genes") %in% colnames(df)))
        df <- df[, c("from", "to", "reversible", "genes")]
    else
        warning("Column names do not match. The first four columns are used for from, to, reversible and genes!")
    vxx <- setdiff(c(df$from, df$to), vnames(g))
    if(! is.empty(vxx)) g <- add.vertices(g, length(vxx), name=vxx)

    robj <- Reactions(g, list.only=FALSE)
    for(i in 1:nrow(df)) {
        ss <- strsplit(df[i, 1], "[ ;,]+")[[1]]
        pp <- strsplit(df[i, 2], "[ ;,]+")[[1]]
        genes <- strsplit(df[i, 4], "[ ;,]+")[[1]]
        ne1 <- length(robj)
        robj <- rset_append(robj, ss, pp, genes)
        ne2 <- length(robj)
        ## add edge if only new reaction added
        if(ne2 > ne1) g <- add.edges(g, expand.grid(ss, pp), reaction=paste0("RX", ne2))
        if(df[i, 3]) {
            robj <- rset_append(robj, pp, ss, genes)
            ne3 <- length(robj)
            if(ne3 > ne2) g <- add.edges(g, expand.grid(pp, ss), reaction=paste0("RX", ne3))
        }
    }
    attr(g, "reactions") <- robj
    g
}

