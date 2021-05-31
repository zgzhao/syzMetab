
#' @title general helper function
#' @description general helpers for retrieving data from keggPATH, reactions list or mgraph object
#' @name helpers_general
#' @aliases Organism Species Compounds Chemicals Reactions Genes
#' @details All functions below are generic functions. Function names are self-explanatory:
#' - pathInfo(object): pathway info
#' - Organism(object): organism name
#' - Species(object): alias of Organism
#' - Genes(object): all genes/orthologs
#' - Compounds(object): all compounds
#' - Chemicals(object): alias of Compounds
#' - Reactions(object, list.only=TRUE): reactions list or ReactionSet object (list.only=FALSE for mgraph object)
#' @examples
#' library(gmetab)
#' d.path <- file.path(path.package("gmetab"), "KEGG")
#' pp <- make_mpath("ko00010", d.path)
#' ##
#' ## KEGG pathway data
#' pathInfo(pp)
#' Organism(pp)
#' Compounds(pp)
#' Chemicals(pp)
#' Genes(pp)
#' # Reactions(pp)
#' ##
#' ## metabolic graph
#' gg <- mgraph_from_mpath(pp)
#' pathInfo(gg)
#' Organism(gg)
#' Compounds(gg)
#' Chemicals(gg)
#' Genes(gg)
#' # Reactions(gg)
#' Reactions(gg, list.only=FALSE)
#' @author ZG Zhao
NULL

#' @export
setGeneric("pathInfo", function(object) standardGeneric("pathInfo"))
#' @export
setGeneric("Genes", function(object) standardGeneric("Genes"))
#' @export
setGeneric("Compounds", function(object) standardGeneric("Compounds"))
Chemicals <- function(...) Compounds(...)
#' @export
setGeneric("Reactions", function(object, ...) standardGeneric("Reactions"))
#' @export
setGeneric("Organism", function(object) standardGeneric("Organism"))
Species <- function(...) Species(...)
#' @export
setGeneric("Substrates", function(object) standardGeneric("Substrates"))
#' @export
setGeneric("Substrates<-", function(object, value) standardGeneric("Substrates<-"))
#' @export
setGeneric("Products", function(object) standardGeneric("Products"))
#' @export
setGeneric("Products<-", function(object, value) standardGeneric("Products<-"))

## ================================================================================
setMethod("pathInfo", "keggPATH", function(object){
    object@pathInfo
})

setMethod("Organism", "keggPATH", function(object){
    object@pathInfo$org
})
setMethod("Organism", "ReactionSet", function(object){
    object@organism
})
setMethod("Organism", "mgraph", function(object){
    x <- attr(object, "reactions")
    x@organism
})

setMethod("Compounds", "keggPATH", function(object){
    rx <- sapply(object@reactions, FUN=function(x) x[c("substrate", "product")])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("Compounds", "ReactionList", function(object){
    rx <- sapply(object, FUN=function(x) x[c("substrate", "product")])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("Compounds", "ReactionSet", function(object){
    rx <- sapply(object@reaction, FUN=function(x) x[c("substrate", "product")])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("Compounds", "mgraph", function(object){
    vnames(object)
})
setMethod("Compounds", "igraph", function(object){
    vnames(object)
})

setMethod("Genes", "keggPATH", function(object){
    rx <- sapply(object@reactions, FUN=function(x) x[["gene"]])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("Genes", "ReactionList", function(object){
    rx <- sapply(object, FUN=function(x) x[["gene"]])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("Genes", "ReactionSet", function(object){
    rx <- sapply(object@reaction, FUN=function(rr) rr[["gene"]])
    rx <- sort(unique(unlist(rx)))
    names(rx) <- NULL
    return(rx)
})
setMethod("Genes", "mgraph", function(object){
    rtns <- Reactions(object)
    rx <- sapply(rtns, FUN=function(rr) rr[["gene"]])
    rx <- unlist(rx)
    sort(unique(rx))
})

setMethod("Reactions", "keggPATH", function(object){
    object@reactions
})
setMethod("Reactions", "ReactionSet", function(object){
    object@reaction
})
setMethod("Reactions", "mgraph", function(object, list.only=TRUE){
    x <- attr(object, "reactions")
    if(list.only) return(x@reaction)
    else return(x)
})

setMethod("Substrates", "mgraph", function(object){
    attr(object, "substrates")
})

setMethod("Products", "mgraph", function(object){
    attr(object, "products")
})

setReplaceMethod("Substrates", c("mgraph", "character"),
                 function(object, value){
    attr(object, "substrates") <- value
    object
})
setReplaceMethod("Products", c("mgraph", "character"),
                 function(object, value){
    attr(object, "products") <- value
    object
})
