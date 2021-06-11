#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: mgraph.R
# Description: handling metabolic graph (mgraph) objects
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-05-31 12:10:39

#' @name make_mgraph
#' @title make graph from various data
#' @aliases make_bgraph make_rgraph make_ggraph
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
    rsobj <- make_rset(object)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "character", function(object, d.path="KEGG"){
    rsobj <- rset_from_kos(object, d.path)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "ReactionList", function(object, org){
    rsobj <- make_rset(object, org)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "list", function(object, org){
    rvalid <- sapply(object, FUN=function(x){
        all(c("substrate", "product", "gene", "reversible") %in% names(x))
    })
    if(! all(rvalid)) stop("Invalid reaction list.")
    rsobj <- make_rset(rtns, org)
    make_mgraph(rsobj)
})

