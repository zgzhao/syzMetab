
#' @title general helper function
#' @description general helpers for retrieving data from KEGGmeta, reactions list or mgraph object
#' @name helpers_general
#' @aliases getOrganism getCPDs getGenes getReactions getAliases Reactions Compounds Species Organism Genes Aliases
#' @details All functions below are generic functions. Function names are self-explanatory:
#' - pathInfo(object)
#' - getOrganism(object)
#' - getGenes(object): all genes/orthologs
#' - getCPDs(object): all compounds
#' - getReactions(object, list.only=TRUE): reactions list or ReactionSet object (list.only=FALSE for mgraph object)
#' - getAliases(object)
#'
#' Alias functions:
#' - Organism(object)
#' - Species(object)
#' - Compounds(object)
#' - Chemicals(object)
#' - Genes(object)
#' - Reactions(object, list.only=TRUE)
#' - Aliases(object)
#' @examples
#' library(gmetab)
#' d.path <- file.path(path.package("gmetab"), "KEGG")
#' kinfo <- kmeta_from_ko("ko00010", d.path)
#' gg <- mgraph_from_kmeta(kinfo)
#' ##
#' pathInfo(kinfo)
#' getOrganism(kinfo)
#' Organism(kinfo)
#' getCPDs(kinfo)
#' Compounds(kinfo)
#' Chemicals(kinfo)
#' getGenes(kinfo)
#' # str(getReactions(kinfo))
#' # str(Reactions(kinfo))
#' ##
#' getOrganism(gg)
#' Organism(gg)
#' getCPDs(gg)
#' Compounds(gg)
#' Chemicals(gg)
#' Genes(gg)
#' # getReactions(gg)
#' # Reactions(gg)
#' Reactions(gg, list.only=FALSE)
#' getAliases(gg)
#' Aliases(gg)
#' @author ZG Zhao
NULL

#' @export
setGeneric("pathInfo", function(object, ...) standardGeneric("pathInfo"))
#' @export
setGeneric("getGenes", function(object, ...) standardGeneric("getGenes"))
#' @export
setGeneric("getCPDs", function(object, ...) standardGeneric("getCPDs"))
#' @export
setGeneric("getReactions", function(object, ...) standardGeneric("getReactions"))
#' @export
setGeneric("getOrganism", function(object) standardGeneric("getOrganism"))
#' @export
setGeneric("getAliases", function(object, ...) standardGeneric("getAliases"))

## aliases
#' @export
Reactions <- function(object, ...) getReactions(object, ...)
#' @export
Aliases <- function(object, ...) getAliases(object, ...)
#' @export
Genes <- function(object, ...) getGenes(object, ...)
#' @export
Compounds <- function(object, ...) getCPDs(object, ...)
#' @export
Chemicals <- function(object, ...) getCPDs(object, ...)
#' @export
Species <- function(object, ...) getOrganism(object, ...)
#' @export
Organism <- function(object, ...) getOrganism(object, ...)

## ================================================================================
setMethod("pathInfo", "KEGGmeta", function(object){
    object@pathInfo
})

setMethod("getOrganism", "KEGGmeta", function(object){
    object@pathInfo$org
})
setMethod("getOrganism", "ReactionSet", function(object){
    object@organism
})
setMethod("getOrganism", "mgraph", function(object){
    x <- attr(object, "reactions")
    x@organism
})

setMethod("getCPDs", "KEGGmeta", function(object){
    rx <- sapply(object@reactions, FUN=function(x) x[c("substrate", "product")])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("getCPDs", "ReactionList", function(object){
    rx <- sapply(object, FUN=function(x) x[c("substrate", "product")])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("getCPDs", "ReactionSet", function(object){
    rx <- sapply(object@reaction, FUN=function(x) x[c("substrate", "product")])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("getCPDs", "mgraph", function(object){
    vnames(object)
})
setMethod("getCPDs", "igraph", function(object){
    vnames(object)
})

setMethod("getGenes", "KEGGmeta", function(object){
    rx <- sapply(object@reactions, FUN=function(x) x[["gene"]])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("getGenes", "ReactionList", function(object){
    rx <- sapply(object, FUN=function(x) x[["gene"]])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("getGenes", "ReactionSet", function(object){
    rx <- sapply(object@reaction, FUN=function(rr) rr[["gene"]])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("getGenes", "mgraph", function(object){
    rtns <- getReactions(object)
    rx <- sapply(rtns, FUN=function(rr) rr[["gene"]])
    rx <- unlist(rx)
    sort(unique(rx))
})

setMethod("getReactions", "KEGGmeta", function(object){
    object@reactions
})
setMethod("getReactions", "ReactionSet", function(object){
    object@reaction
})
setMethod("getReactions", "mgraph", function(object, list.only=TRUE){
    x <- attr(object, "reactions")
    if(list.only) return(x@reaction)
    else return(x)
})

setMethod("getAliases", "ReactionSet", function(object){
    rx <- object@alias
    rx
})
setMethod("getAliases", "mgraph", function(object){
    x <- attr(object, "reactions")
    rx <- x@alias
    rx
})
