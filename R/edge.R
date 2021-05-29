#' get edge names
#'
#' Similar functions: \code{\link{vnames}}, \code{\link{enames}}, \code{\link{rnames}}, \code{\link{vcount}}, \code{\link{ecount}}
#' @title edge names
#' @param object mgraph object
#' @author ZG Zhao
#' @export
setGeneric("enames", function(object) standardGeneric("enames"))
setMethod("enames", "mgraph", function(object){
    as_ids(E(object))
})


#' Get or set edge attribute/data
#'
#' Friendly function for edge data/attributes retrieving or setting. See "Example".
#' @title edge data
#' @param g igraph/mgraph object
#' @param a.name character, attribute name
#' @param e.names vector of character (edge names) or integer (indices)
#' @seealso \code{\link{vdata}}
#' @author ZG Zhao
#' @examples
#' library(mgraph)
#' d.path <- file.path(path.package("mgraph"), "KEGG")
#' gg <- mgraph_from_kos("ko00010", d.path)
#' ## igraph style
#' E(gg)$reaction
#' E(gg)$reaction[1:3]
#' E(gg)$MP <- 0
#' ## mgraph style
#' edata(gg, "reaction")
#' edata(gg, "reaction", 1:3)
#' (ee <- enames(gg))
#' edata(gg, "reaction", ee[1:3])
#' edata(gg, "MP") <- 1
#' edata(gg, "MP")
#' @export
edata <- function(g, a.name, e.names) {
    if(! is.xgraph(g)) return(NULL)
    a.name <- unlist(a.name)[1]
    ## data presented in MRreactions
    xnames <- c("substrate", "product", "gene")
    if(is.mgraph(g) && a.name %in% xnames) {
        rtns <- Reactions(g)
        rx <- sapply(rtns, FUN=function(rr) rr[[a.name]])
        if(missing(e.names)) return(rx)

        ss <- (1:ecount(g) %in% e.names) | (enames(g) %in% e.names)
        r.names <- E(g)$reaction
        r.names <- r.names[ss]
        if(sum(ss) == 1) return(rx[[r.names]])
        else return(rx[r.names])
    }

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
        stop("Cannot modify this attribute!")
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


#' @export
outEdges <- function(es, g) adj(g, es)

#' @export
edgesNames <- function(g){
    names(edgeData(g))
}


#' @export
edgesInfo <- function(g, attr.name, edge.names = NULL) {
    ess <- edgesNames(g)
    if(is.empty(ess)) return(NULL)

    edge.names <- intersect(ess, edge.names)
    if(length(edge.names) < 1) edge.names <- ess

    res <- lapply(edge.names, FUN = function(x){
        ee <- strsplit(x, "|", fix=TRUE)[[1]]
        xx <- unlist(edgeData(g, ee[1], ee[2], attr.name))
        names(xx) <- NULL
        xx
    })
    if(all(lapply(res, length) == 1)) {
        res <- unlist(res)
    }
    names(res) <- edge.names
    res
}

#' @export
deleteEdges <- function(g, ess) {
    if(is.matrix(ess) || is.data.frame(ess)) {
        ss <- unlist(ess[, 1])
        ee <- unlist(ess[, 2])
    } else {
        ess <- intersect(ess, edgesNames(g))
        if(length(ess) < 1) return(g)
        ess <- sapply(ess, FUN = function(x) strsplit(x, "|", fix=TRUE)[[1]])
        ss <- unlist(ess[1, ])
        ee <- unlist(ess[2, ])
    }
    for(i in 1:length(ss)) {
        if(is.adjacentVs(g, ss[i], ee[i], "out"))
            g <- removeEdge(ss[i], ee[i], g)
    }
    g
}

.edgesInPaths <- function(spp, from, to){
    ee <- lapply(spp, FUN = function(x) {
        aa <- x %in% from
        bb <- x %in% to
        if(sum(aa, bb) < 2) return(NULL)
        ## FIXME: more than one source/target in path
        aa <- which(aa)
        bb <- which(bb)
        if(aa > bb) return(NULL)
        x <- x[aa:bb]
        n <- length(x)
        paste(x[1:(n - 1)], x[2:n], sep="|")
    })

    unique(unlist(ee))
}
