#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: node.R
# Description: handling node/vertex associates
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2021-06-02 08:37:18

#' @export
setGeneric("vnames", function(object) standardGeneric("vnames"))
setMethod("vnames", "igraph", function(object){
    as_ids(V(object))
})
setMethod("vnames", "xgraph", function(object){
    as_ids(V(object))
})

#' @export
vdata <- function(g, a.name, v.names) {
    if(! is.xgraph(g)) return(NULL)
    a.name <- a.name[1]
    ee <- substitute(V(g)$`x`, list(x=a.name))
    rx <- eval(ee)
    if(is.null(rx) || missing(v.names)) return(rx)
    ss <- (1:vcount(g) %in% v.names) | (vnames(g) %in% v.names)
    rx[v.names]
}

#' @export
`vdata<-` <- function(g, a.name, v.names, value) {
    a.name <- a.name[1]
    if(missing(v.names)) {
        ee <- substitute(V(g)$`attx` <- value, list(attx=a.name, value=value))
    } else {
        ss <- (1:vcount(g) %in% v.names) | (vnames(g) %in% v.names)
        ndx <- which(ss)
        ee <- substitute(V(g)$`attx`[n] <- value, list(attx=a.name, n=ndx, value=value))
    }
    eval(ee)
    g
}

#' delete vertices/nodes from igraph or mgraph object
#'
#' Same as functions in igraph package exception for retaining graph attributes. Refer to `?igraph::delete.vertices`, `?igraph::delete.vertices` for details.
#' @title delete vertices
#' @aliases delete_vertices add.vertices add_vertices
#' @return igraph/mgraph object
#' @author ZG Zhao
#' @export
setGeneric("delete.vertices", function(object, vs) standardGeneric("delete.vertices"))
#' @export
delete_vertices <- function(...) delete.vertices(...)
setMethod("delete.vertices", "igraph", function(object, vs) {
    g <- igraph::delete.vertices(object, vs)
    attributes(g) <- attributes(object)
    g
})

setMethod("delete.vertices", "xgraph", function(object, vs) {
    g <- igraph::delete.vertices(object, vs)
    attributes(g) <- attributes(object)
    g
})

#' @export
setGeneric("add.vertices", function(object, x, ...) standardGeneric("add.vertices"))
#' @export
add_vertices <- function(...) add.vertices(...)
setMethod("add.vertices", c("igraph", "numeric"), function(object, x, ...) {
    g <- igraph::add.vertices(object, x, ...)
    attributes(g) <- attributes(object)
    g
})
setMethod("add.vertices", c("igraph", "character"), function(object, x, ...) {
    g <- igraph::add.vertices(object, length(x), name=x, ...)
    attributes(g) <- attributes(object)
    g
})
setMethod("add.vertices", "xgraph", function(object, x, ...) {
    g <- add.vertices(as(object, "igraph"), x, ...)
    attributes(g) <- attributes(object)
    g
})


#' Get all vertices assessed by given vertices.
#'
#' "Accessed" means that not all the distances from (mode="out") or to (mode="in") given vertices are infinite.
#' @title vertices assessed by given vertices
#' @param g graph object
#' @param v.names character vector, names of vertices
#' @param mode character, "out" (default), "in" or "all", direction from given vertices.
#' @return vertices vector
#' @author ZG Zhao
#' @export
vs_accessed_by <- function(g, v.names, mode=c("out", "in", "all")){
    vns <- vnames(g)
    v.names <- intersect(v.names, vns)
    if(is.empty(v.names)) return(NULL)
    dd <- distances(g, which(vns %in% v.names), mode=mode[1])
    ss <- apply(dd, 2, FUN=function(x) ! all(is.infinite(x)))
    names(ss)[ss]
}

#' @export
vs_adjacent <- function(g, v.names, mode = c("out", "out", "both")) {
    v.names <- intersect(v.names, vnames(g))
    lapply(adjacent_vertices(g, v.names, mode=mode[1]), names)
}

#' test whether nodes are connected (have any route)
#'
#' Note that non-connecting nodes are judged by infinite distance.
#' @title test connected
#' @param g graph object
#' @param s source node(s)
#' @param p target edge(s)
#' @return logi, TRUE if any target is accessed by any source
#' @author ZG Zhao
#' @export
vis_connected <- function(g, s, p) {
    if(vcount(g) < 1) return(FALSE)
    if(is.empty(s) || is.empty(p)) return(FALSE)
    any(p %in% vs_accessed_by(g, s, "out"))
}
