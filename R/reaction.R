#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: reaction.R
# Description: handling reaction associates
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-06-02 08:31:27

#' @name rnames
#' @title names of graph entities
#' @description Get names (instead of ids) of nodes, edges or reactions (xgraph only).
#' @details
#' - enames: edges names (igraph)
#' - vnames: vertex/node names (igraph)
#' - rnames: For graph object, it returns the reaction names (may be duplicated) along vnames. For "ReactionSet" object, they are the (unique) names of reactions.
#' @aliases vnames enames
#' @param g graph (igraph or xgraph) object
#' @param xg xgraph object
#' @param rset ReactionSet object
#' @usage
#' enames(g)
#' vnames(g)
#' rnames(xg)
#' rnames(rset)
#' @return character vector
#' @author ZG Zhao
NULL

#' @name rcount
#' @title count of graph entities
#' @description Get count/number of nodes, edges or reactions
#' @details
#' - ecount: edges names (igraph)
#' - vcount: vertex/node names (igraph)
#' - rcount: number of unique reactions
#' @aliases vcount ecount
#' @param g graph (igraph or xgraph) object
#' @param xg xgraph object
#' @param rset ReactionSet object
#' @usage
#' ecount(g)
#' vcount(g)
#' rcount(xg)
#' rcount(rset)
#' @return character vector
#' @author ZG Zhao
NULL

#' @name rdata
#' @title get/set attribute data
#' @aliases edata vdata
#' @description Helper functions attribute data
#' @details Provide alternative ways to get or set graph attributes.
#' @param g igraph/mgraph object
#' @param a.name attribute name
#' @param v.names (optional) node names or indices
#' @param e.names (optional) edge names or indices
#' @param r.names (optional) reaction names (for subseting)
#' @usage
#' vdata(g, a.name, v.name)
#' edata(g, a.name, e.name)
#' rdata(g, a.name, r.name)
#' vdata(g, a.name) <- vals
#' edata(g, a.name) <- vals
#' rdata(g, a.name) <- vals
#' @author ZG Zhao
NULL


#' @export
setGeneric("rnames", function(object) standardGeneric("rnames"))
setMethod("rnames", "xgraph", function(object){
    E(object)$reaction
})
setMethod("rnames", "ReactionSet", function(object){
    names(Reactions(object))
})

#' @export
setGeneric("rcount", function(object) standardGeneric("rcount"))
setMethod("rcount", "xgraph", function(object){
    length(Reactions(object))
})
setMethod("rcount", "ReactionSet", function(object){
    length(Reactions(object))
})

#' @export
rdata <- function(g, a.name, x.names) {
    if(! is.mgraph(g)) return(NULL)
    a.name <- unlist(a.name)[1]
    if(tolower(a.name) == "organism")
        return(Organism(g))

    rtns <- Reactions(g)
    rx <- lapply(rtns, FUN=function(rr) rr[[a.name]])
    if(missing(x.names)) return(rx)

    ss <- (1:length(rx) %in% x.names) | (names(rx) %in% x.names)
    x.names <- names(rx)[ss]
    if(sum(ss) < 1) return(NULL)
    else if(sum(ss) == 1) return(rx[[x.names]])
    else return(rx[x.names])
}

## internal use only!!
`rdata<-` <- function(g, a.name, value) {
    a.name <- a.name[1]
    if(! a.name %in% c("reaction", "organism", "alias"))
        stop("Attributes must be: reaction, organism or alias.")
    rtns <- Reactions(g, list.only=FALSE)
    ee <- substitute(rtns@`attx` <- value, list(attx=a.name, value=value))
    eval(ee)
    attr(g, "reactions") <- rtns
    g
}


#' Mostly for internal use.
#'
#' "ReactionSet" is S4 class designed for mgraph. Reactions are presented in different forms from those of \code{\link{keggPATH}}.
#' @title generate ReactionSet object
#' @param object keggPATH or list object. If feeding a list object, make sure each list elements are also lists containing substrate, product, gene and reversible elements.
#' @param org character, "ko" or organism identifier.
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @param ... other params
#' @return ReactionSet object
#' @author ZG Zhao
#' @export
setGeneric("make_rset", function(object, org="ko", d.path="KEGG", ...) standardGeneric("make_rset"))
setMethod("make_rset", "ReactionList", function(object, org){
    rtns <- list()
    hasRev <- "reversible" %in% names(object[[1]])
    for(rx in object) {
        ss <- rx[["substrate"]]
        pp <- rx[["product"]]
        genes <- rx[["gene"]]
        ids <- rx[["name"]]
        revx <- if(hasRev) rx[["reversible"]] else FALSE
        sspp <- expand.grid(ss, pp)
        sspp[[1]] <- as.character(sspp[[1]])
        sspp[[2]] <- as.character(sspp[[2]])
        for(i in 1:nrow(sspp)) {
            ss <- sspp[[1]][i]
            pp <- sspp[[2]][i]
            rtns <- rx_append(rtns, ss, pp, genes, ids)
            if(revx) rtns <- rx_append(rtns, pp, ss, genes, ids)
        }
    }
    names(rtns) <- paste0("RX", 1:length(rtns))
    new("ReactionSet", rtns, org)
})
setMethod("make_rset", "list", function(object, org){
    class(object) <- "ReactionList"
    make_rset(object)
})

setMethod("make_rset", "keggPATH", function(object){
    org <- pathInfo(object)$org
    rlist <- Reactions(object)
    make_rset(rlist, org)
})
## from kos
setMethod("make_rset", "character", function(object, d.path){
    org <- unique(gsub("[0-9]+", "", object))
    org <- setdiff(org, c("", "ko"))
    if(length(org) > 1) stop("Only one organism is allowed.")
    if(length(org) < 1) org <- "ko"
    kndx <- unique(gsub("[a-z]+", org, object))
    rtns <- list()
    for(kx in kndx) {
        xinfo <- make_mpath(kx, d.path)
        rtns <- c(rtns, Reactions(xinfo))
    }
    make_rset(rtns, org)
})

#' Append reaction to reaction lists
#'
#' Useful for adding to network spontaneous reactions without associated genes.
#' @title append reaction
#' @param object ReactionList or equivalent list
#' @param x character (substrate name) or data.frame with at least these columns: from, to, genes, reversible
#' @param p character, product name (ignore if x is a data.frame)
#' @param genes character vector, names of genes involved in the reaction  (ignore if x is a data.frame)
#' @param ids additional identifers for reaction  (ignore if x is a data.frame)
#' @param reversible TRUE/FALSE (default). If TRUE, add reactions in both directions. Ignore if x is a data.frame.
#' @return  ReactionList object
#' @author ZG Zhao
#' @export
setGeneric("rlist_append", function(object, x, p, genes, ids=NULL, reversible=FALSE) standardGeneric("rlist_append"))
setMethod("rlist_append", c("ReactionList", "character"), function(object, x, p, genes, ids, reversible){
    rlist <- rx_append(object, x, p, genes, ids)
    if(reversible) rlist <- rx_append(rlist, p, x, genes, ids)
    rlist
})
setMethod("rlist_append", c("list", "character"), function(object, x, p, genes, ids, reversible){
    class(object) <- "ReactionList"
    rlist_append(object, x, p, genes, ids, reversible)
})
setMethod("rlist_append", c("ReactionList", "data.frame"), function(object, x){
    for(i in 1:nrow(x))
        object <- rlist_append(object, x$from[i], x$to[i], x$genes[i], reversible = x$reversible[i])
    object
})
setMethod("rlist_append", c("list", "data.frame"), function(object, x){
    class(object) <- "ReactionList"
    rlist_append(object, x)
})

rx_append <- function(rlist, s, p, genes, ids=NULL) {
    s <- sort(s)
    p <- sort(p)
    rx <- list(substrate=s, product=p, gene=genes, ids=ids)

    if(length(rlist) < 1) {
        rlist[["RX1"]] <- rx
    } else {
        r.names <- sapply(rlist, FUN=function(rx){
            ss <- sort(rx[["substrate"]])
            pp <- sort(rx[["product"]])
            paste(c(ss, pp), collapse=" ")
        })
        names(r.names) <- NULL
        xname <- paste(c(s, p), collapse=" ")
        if(xname %in% r.names) {
            ndx <- which(r.names == xname)
            genes <- c(genes, rlist[[ndx]][["gene"]])
            rlist[[ndx]][["gene"]] <- unique(genes)
            ids <- c(ids, rlist[[ndx]][["ids"]])
            rlist[[ndx]][["ids"]] <- unique(ids)
        } else {
            xndx <- paste0("RX", length(rlist) + 1)
            rlist[[xndx]] <- rx
        }
    }
    class(rlist) <- "ReactionList"
    rlist
}
