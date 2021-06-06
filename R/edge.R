#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: edge.R
# Description: handling edge associates
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-06-02 08:39:08


#' @export
setGeneric("enames", function(object) standardGeneric("enames"))
setMethod("enames", "igraph", function(object){
    as_ids(E(object))
})
setMethod("enames", "xgraph", function(object){
    as_ids(E(object))
})


#' @export
edata <- function(g, a.name, e.names) {
    if(! is.xgraph(g)) return(NULL)
    a.name <- unlist(a.name)[1]
    if(a.name %in% c("substrate", "product", "gene"))
        stop("Please use `rdata` for reaction-related data.")
    ## data presented as graph attribute
    ee <- substitute(E(g)$`x`, list(x=a.name))
    rx <- eval(ee)
    if(is.null(rx) || missing(e.names)) return(rx)
    ss <- (1:ecount(g) %in% e.names) | (enames(g) %in% e.names)
    rx[ss]
}

#' @export
`edata<-` <- function(g, a.name, e.names, value) {
    a.name <- a.name[1]
    if(a.name %in% c("substrate", "product", "gene"))
        stop("Please use `rdata` for reaction-related data.")
    if(missing(e.names)) {
        ee <- substitute(E(g)$`attx` <- value, list(attx=a.name, value=value))
    } else {
        ss <- (1:ecount(g) %in% e.names) | (enames(g) %in% e.names)
        ndx <- which(ss)
        ee <- substitute(E(g)$`attx`[n] <- value, list(attx=a.name, n=ndx, value=value))
    }
    eval(ee)
    g
}


#' delete edges from igraph or mgraph object
#'
#' Same as functions in igraph package exception for retaining graph attributes. Refer to `?igraph::delete.edges`
#' @title delete edges
#' @aliases delete_edges
#' @param object igraph/mgraph object
#' @param es vector: edge ids (integer) or names (character)
#' @return igraph/mgraph object
#' @author ZG Zhao
#' @export
setGeneric("delete.edges", function(object, es) standardGeneric("delete.edges"))
#' @export
delete_edges <- function(...) delete.edges(...)

setMethod("delete.edges", "igraph", function(object, es) {
    g <- igraph::delete.edges(object, es)
    attributes(g) <- attributes(object)
    g
})
setMethod("delete.edges", "xgraph", function(object, es) {
    ss1 <- 1:ecount(object) %in% es
    ss2 <- enames(object) %in% es
    ss3 <- rnames(object) %in% es
    g <- igraph::delete.edges(object, which(ss1 | ss2 | ss3))
    attributes(g) <- attributes(object)
    g
})

#' add edges for igraph or mgraph object
#'
#' Versatile add edges fucntion: accept sequences (like igraph::add.edges), vertex names or data.frame (first 2 columns are `from` and `to`). Refer to `igraph::add.edges` for other details.
#' @title add edges
#' @aliases add_edges
#' @param object igraph/mgraph object
#' @param edges various data type
#' - sequence of integer or character: like igraph add.edges
#' - edge names: each name is a pair of nodes seperated by "|" such as "A|B"
#' - data.frame or matrix: first 2 columns are treated as sources and targets.
#' @return graph object
#' @author ZG Zhao
#' @export
setGeneric("add.edges", function(object, edges, ...) standardGeneric("add.edges"))
#' @export
add_edges <- function(...) add.edges(...)

setMethod("add.edges", c("igraph", "vector"), function(object, edges, ...) {
    if(all(grepl("|", edges, fixed=TRUE))) {
        ess <- sapply(edges, FUN=function(x) strsplit(x, "|", fixed=T)[[1]])
        edges <- as.matrix(ess)
    }
    ss <- !(edges %in% 1:vcount(object) | edges %in% vnames(object))
    if(sum(ss) > 0) object <- add.vertices(object, edges[ss])
    igraph::add.edges(object, edges, ...)
})
setMethod("add.edges", c("igraph", "matrix"), function(object, edges, ...) {
    vss <- setdiff(c(edges[, 1], edges[, 2]), vnames(object))
    if(!is.empty(vss)) object <- add.vertices(object, vss)

    ess <- apply(edges[, 1:2], 1, paste, collapse="|")
    tt <- ! ess %in% enames(object)
    if(sum(tt) < 1) return(object)
    ess <- ess[tt]
    ## call add.edges(object, vector, ...)
    add.edges(object, ess, ...)
})
setMethod("add.edges", c("igraph", "data.frame"), function(object, edges, ...) {
    vss <- setdiff(unlist(edges), vnames(object))
    if(!is.empty(vss)) object <- add.vertices(object, vss)

    ess <- apply(edges[, 1:2], 1, paste, collapse="|")
    tt <- ! ess %in% enames(object)
    if(sum(tt) < 1) return(object)
    ess <- ess[tt]
    ## call add.edges(object, vector, ...)
    add.edges(object, ess, ...)
})
setMethod("add.edges", "xgraph", function(object, edges, ...) {
    g <- as(object, "igraph")
    g <- add.edges(g, edges, ...)
    attributes(g) <- attributes(object)
    g
})

