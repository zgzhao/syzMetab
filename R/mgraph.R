#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: mgraph.R
# Description: handling metabolic graph (mgraph) objects
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-05-31 12:10:39

#' @name make_mgraph
#' @title make graph from various data
#' @aliases make_bgraph make_rgraph make_ggraph make_bgraph
#' @description Make metabolic graph from various type of data: KDataSet, ReactionSet, ReactionList, list, KOs or other xgraph object.
#' @details Metabolic network can be represented by metabolite, reaction, gene graph or st-cuts of gene sets.
#' - make_mgraph (make metabolite graph): nodes are metabolites and edges are reactions identified by gene sets
#' - make_rgraph (make reaction graph): nodes are reactions
#' - make_ggraph (make gene graph): nodes are genes
#' - make_bgraph (make bipartite graph): derived from st-cuts of mgraph. Nodes are genes and st-cuts of gene sets. Edges are links between gene sets and the genes included.
#' @usage
#' - make_mgraph(object, ...): from KDataSet, ResctionSet, ReactionList, list or character vector (kos)
#' - make_rgraph(object, ...): from KDataSet, ResctionSet, ReactionList or character vector (kos)
#' - make_ggraph(object, ...): from KDataSet, ResctionSet, ReactionList or character vector (kos)
#' - make_bgraph(object, ...): from mgraph or stcuts (refer to \code{\link{make_stcuts}})
#' @param object KDataSet, ReactionSet ReactionList, character vector or list.
#' @param ... additional parameters depends on "object"
#' - org: organism/species indentifier, parameter for "object" without organism info such as ReactionList and list
#' - d.path: character, refer to \code{\link{KEGG_get}} for detail. Used when "object" is a character vector (kos).
#' @author zhao
#' @examples
#' ## NOT RUN
#' library(syzMetab)
#' pp <- make_kdset("ko00010", d.path="KEGG")
#' make_mgraph(pp)
#' make_mgraph(Reactions(pp), org="ko")
#' make_mgraph("ko00010")
#' gm <- make_mgraph(c("ko00010", "ko00020"))
#' plot(gm)
#' gr <- make_rgraph("ko00010")
#' plot(gr)
#' gg <- make_ggraph("ko00010")
#' plot(gg)
#' chem1 <- c("C00031", "C00221", "C00267", "C01172", "C01451", "C06186")
#' chem2 <- "C00022"
#' gb <- make_bgraph(gm, s=chem1, t=chem2)
#' plot(gb)
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
setMethod("make_mgraph", "KDataSet", function(object){
    rsobj <- make_rset(object)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "character", function(object, d.path="KEGG"){
    rsobj <- make_rset(object, d.path)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "ReactionList", function(object, org="ko"){
    rsobj <- make_rset(object, org)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "list", function(object, org="ko"){
    rvalid <- sapply(object, FUN=function(x){
        all(c("substrate", "product", "gene", "reversible") %in% names(x))
    })
    if(! all(rvalid)) stop("Invalid reaction list.")
    rsobj <- make_rset(object, org)
    make_mgraph(rsobj)
})

