#' get reaction names
#'
#' Similar functions: \code{\link{vnames}}, \code{\link{enames}}, \code{\link{rnames}}, \code{\link{vcount}}, \code{\link{ecount}}
#' @title reaction names
#' @param object mgraph or RSet object
#' @author ZG Zhao
#' @export
setGeneric("rnames", function(object) standardGeneric("rnames"))
setMethod("rnames", "mgraph", function(object){
    E(object)$reaction
})
setMethod("rnames", "RSet", function(object){
    names(Reactions(object))
})

#' @export
setGeneric("rcount", function(object) standardGeneric("rcount"))
setMethod("rcount", "mgraph", function(object){
    length(Reactions(object))
})
setMethod("rcount", "RSet", function(object){
    length(Reactions(object))
})

#' Get or set reaction data
#'
#' Friendly function for edge data/attributes retrieving or setting. See "Example".
#' @title reaction data
#' @param g mgraph or RSet object
#' @param a.name character, attribute name
#' @param x.names vector of character (reaction names) or integer (ids)
#' @seealso \code{\link{vdata}}
#' @author ZG Zhao
#' @export
rdata <- function(g, a.name, x.names) {
    if(! is.mgraph(g)) return(NULL)
    a.name <- unlist(a.name)[1]
    if(tolower(a.name) == "organism") Organism(g)
    else if(tolower(a.name) == "alias") Aliases(g)

    rtns <- Reactions(g)
    rx <- lapply(rtns, FUN=function(rr) rr[[a.name]])
    if(missing(x.names)) return(rx)

    ss <- (1:length(rx) %in% x.names) | (names(rx) %in% x.names)
    x.names <- names(rx)[ss]
    if(sum(ss) < 1) return(NULL)
    else if(sum(ss) == 1) return(rx[[x.names]])
    else return(rx[x.names])
}

#' @export
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
#' "RSet" is S4 class designed for mgraph. Reactions are presented in different forms from those of \code{\link{KEGGmeta}}.
#' @title generate RSet object
#' @param kinfo KEGGmeta or list object. If feeding a list object, make sure each list elements are also lists containing substrate, product, gene and reversible elements.
#' @param org character, a KEGG organism string. Use if only kinfo is a plain list.
#' @return RSet object
#' @author ZG Zhao
#' @export
as_rset <- function(kinfo, org="ko") {
    if(is.kmeta(kinfo)) {
        org <- pathInfo(kinfo)$org
        kinfo <- getReactions(kinfo)
    }
    rtns <- list()
    for(rx in kinfo) {
        ss <- sort(rx[["substrate"]])
        pp <- sort(rx[["product"]])
        genes <- rx[["gene"]]
        rtns <- rlist_append(rtns, ss, pp, genes)
        if(rx[["reversible"]])
            rtns <- rlist_append(rtns, pp, ss, genes)
    }
    names(rtns) <- paste0("RX", 1:length(rtns))
    new("RSet", rtns, org)
}

#' Append reaction to RSet object
#'
#'
#' @title append reaction
#' @param robj RSet object
#' @param s character, name of substrate
#' @param p character, name of product
#' @param genes character vector, names of genes involved in the reaction
#' @return  RSet object
#' @author ZG Zhao
#' @export
rset_append <- function(robj, s, p, genes) {
    rlist <- rlist_append(robj@reaction, s, p, genes)
    robj@reaction <- rlist
    robj
}

#' @export
rset_merge_chems <- function(robj, chems, alias) {
    chems <- intersect(Chemicals(robj), chems)
    if(is.empty(chems)) return(robj)

    alias <- alias[1]
    xalias <- Aliases(robj)
    xalias[alias] <- chems
    rlist <- Reactions(robj)
    rlist <- lapply(rlist, FUN=function(rx){
        chem1 <- rx[["substrate"]]
        chem1[chem1 %in% chems] <- alias
        chem1 <- unique(chem1)
        chem2 <- rx[["substrate"]]
        chem2[chem2 %in% chems] <- alias
        chem2 <- unique(chem2)
        if(alias %in% chem1 && alias %in% chem2) {
            chem1 <- setdiff(chem1, alias)
            chem2 <- setdiff(chem2, alias)
        }
        if(is.empty(chem1) || is.empty(chem2)) return(NA)
        rx[["substrate"]] <- chem1
        rx[["product"]] <- chem2
        rx[["reversible"]] <- FALSE
        rx
    })
    ss <- sapply(rlist, FUN=is.empty)
    rlist <- rlist[!ss]
    ## rebuild object: merge some reactions
    robj <- as_rset(rlist, Organism(robj))
    robj@alias <- xalias
    robj
}

rlist_append <- function(rlist, s, p, genes) {
    s <- sort(s)
    p <- sort(p)
    rx <- list(substrate=s, product=p, gene=genes)

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
        } else {
            xndx <- paste0("RX", length(rlist) + 1)
            rlist[[xndx]] <- rx
        }
    }
    rlist
}
