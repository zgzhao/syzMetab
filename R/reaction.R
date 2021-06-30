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
#' "ReactionSet" is S4 class designed for mgraph. Reactions are presented in different forms from those of \code{\link{KDataSet}}.
#' @title generate ReactionSet object
#' @param object KDataSet or list object. If feeding a list object, make sure each list elements are also lists containing substrate, product, gene and reversible elements.
#' @param org character, "ko" or organism identifier.
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @param ... other params
#' @return ReactionSet object
#' @author ZG Zhao
#' @export
setGeneric("make_rset", function(object, org="ko", d.path="KEGG", ...) standardGeneric("make_rset"))
setMethod("make_rset", "list", function(object, org){
    rtns <- list()
    for(rx in object) {
        aa <- rx[["substrate"]]
        bb <- rx[["product"]]
        cc <- rx[["gene"]]
        dd <- if("name" %in% names(rx)) rx[["name"]] else NULL
        ee <- if("reversible" %in% names(rx)) rx[["reversible"]] else FALSE
        rtns <- add.reactions(rtns, ds=aa, t=bb, genes=cc, name=dd, reversible=ee)
    }
    names(rtns) <- paste0("RX", 1:length(rtns))
    new("ReactionSet", rtns, org)
})
setMethod("make_rset", "ReactionList", function(object, org){
    class(object) <- "list"
    make_rset(object, org)
})

setMethod("make_rset", "KDataSet", function(object){
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
        xinfo <- make_kdset(kx, d.path)
        rtns <- c(rtns, Reactions(xinfo))
    }
    make_rset(rtns, org)
})

#' Append reaction to reaction lists
#'
#' Useful for adding to network spontaneous reactions without associated genes.
#' @title append reaction
#' @param object ReactionList or equivalent list
#' @param ds character (substrate name) or data.frame with at least these columns: from, to, genes, reversible
#' @param t character, product name (ignore if ds is a data.frame)
#' @param genes character vector, names of genes involved in the reaction  (ignore if ds is a data.frame)
#' @param name additional identifer(s) for reaction  (ignore if ds is a data.frame)
#' @param reversible TRUE/FALSE (default). If TRUE, add reactions in both directions. Ignore if ds is a data.frame.
#' @return  ReactionList object
#' @author ZG Zhao
#' @export
setGeneric("add.reactions", function(object, ds, t, genes, name=NULL, reversible=FALSE) standardGeneric("add.reactions"))
setMethod("add.reactions", c("list", "character", "character", "character"),
          function(object, ds, t, genes, name, reversible){
    ## All gene names in ReactionSet are normalized.
    genes <- .setGeneNames(unique(genes))
    ## reaction names by s-t names
    r.names <- sapply(object, FUN=function(rx){
        aa <- sort(rx[["substrate"]])
        bb <- sort(rx[["product"]])
        paste(aa, bb, sep=" ")
    })
    names(r.names) <- NULL
    ## take reversibility into account
    sts <- expand.grid(ds, t)
    if(reversible) sts <- rbind(sts, expand.grid(t, ds))
    sts[[1]] <- as.character(sts[[1]])
    sts[[2]] <- as.character(sts[[2]])

    for(ic in 1:nrow(sts)) {
        ss <- sts[ic, 1]
        tt <- sts[ic, 2]
        if(ss == tt) next
        x.name <- paste(ss, tt, sep=" ")

        if(x.name %in% r.names) {
            ndx <- which(r.names == x.name)
            ndx <- ndx[1]  ## one reaction only each time
            genes <- c(genes, object[[ndx]][["gene"]])
            genes <- setdiff(genes, c(NA, ""))
            name <- c(name, object[[ndx]][["name"]])
            object[[ndx]][["gene"]] <- unique(genes)
            object[[ndx]][["name"]] <- unique(name)
        } else {
            rx <- list(substrate=ss, product=tt, gene=genes, name=name)
            xndx <- paste0("RX", length(object) + 1)
            object[[xndx]] <- rx
            r.names <- c(r.names, x.name)
        }
    }
    object
})
setMethod("add.reactions", c("list", "data.frame"), function(object, ds){
    c.names <- colnames(ds)
    if(! all(c("from", "to", "genes", "reversible") %in% c.names))
        stop("Data should contain at least 4 columns: from, to, genes and reversible.")
    for(ic in 1:nrow(ds)) {
        s <- strsplit(ds$from[ic], "[ ;,\\|]+")[[1]]
        t <- strsplit(ds$to[ic], "[ ;,\\|]+")[[1]]
        genes <- strsplit(ds$genes[ic], "[ ;,\\|]+")[[1]]
        xrev <- ds$reversible[ic]
        xid <- if("name" %in% c.names) NULL else ds$name[ic]
        object <- add.reactions(object, ds=s, t=t, genes=genes, name=xid, reversible=xrev)
    }
    object
})
setMethod("add.reactions", c("ReactionSet", "character", "character", "character"),
          function(object, ds, t, genes, name, reversible){
    rlist <- Reactions(object)
    rlist <- add.reactions(rlist, ds, t, genes, name, reversible)
    object@reaction <- rlist
    object
})
setMethod("add.reactions", c("ReactionSet", "data.frame"), function(object, ds){
    rlist <- Reactions(object)
    rlist <- add.reactions(rlist, ds)
    object@reaction <- rlist
    object
})

setMethod("add.reactions", c("mgraph", "character", "character", "character"),
          function(object, ds, t, genes, name, reversible){
    if(Organism(object) != "ko") genes <- .setGeneNames(genes)
    rlist <- Reactions(object, list.only=TRUE)
    rlist <- add.reactions(rlist, ds, t, genes, name, reversible)
    rset <- Reactions(object, list.only=FALSE)
    rset@reaction <- rlist
    ## NOTE lazy method: make mgraph from ReactionSet
    ## TODO: 1. add vertices; 2. add edges with reaction name and gene mappings
    g <- make_mgraph(rset)
    g
})

setMethod("add.reactions", c("mgraph", "data.frame"), function(object, ds){
    c.names <- colnames(ds)
    if(! all(c("from", "to", "genes", "reversible") %in% c.names))
        stop("Data should contain at least 4 columns: from, to, genes and reversible.")
    if(Organism(object) != "ko") ds$genes <- .setGeneNames(ds$genes)
    rlist <- Reactions(object, list.only=TRUE)
    rlist <- add.reactions(rlist, ds)
    rset <- Reactions(object, list.only=FALSE)
    rset@reaction <- rlist
    g <- make_mgraph(rset)
    g
})

#' @export
setGeneric("merge_chems", function(object, chems, nickname) standardGeneric("merge_chems"))
setMethod("merge_chems", "mgraph", function(object, chems, nickname){
    nickname <- nickname[1]
    rlist <- Reactions(object)
    rlist <- lapply(rlist, FUN=function(rx){
        ss <- rx[["substrate"]]
        tt <- rx[["product"]]
        ss[ss %in% chems] <- nickname
        tt[tt %in% chems] <- nickname
        rx[["substrate"]] <- ss
        rx[["product"]] <- tt
        rx
    })
    rtns <- make_rset(rlist, Organism(object))
    make_mgraph(rtns)
})
