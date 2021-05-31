
#' get edge names
#'
#' Similar functions: \code{\link{vnames}}, \code{\link{enames}}, \code{\link{rnames}}, \code{\link{vcount}}, \code{\link{ecount}}
#' @title edge names
#' @param object mgraph object
#' @author ZG Zhao
#' @export
setGeneric("enames", function(object) standardGeneric("enames"))
setMethod("enames", "igraph", function(object){
    as_ids(E(object))
})
setMethod("enames", "xgraph", function(object){
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
#' library(gmetab)
#' d.path <- file.path(path.package("gmetab"), "KEGG")
#' gg <- mgraph_from_kos("ko00010", d.path)
#' ## igraph style
#' E(gg)$reaction
#' E(gg)$reaction[1:3]
#' E(gg)$MP <- 0
#' ## gmetab style
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
#' Refer to \code{\link{igraph::delete.edges}}
#' @title delete edges
#' @aliases delete_edges delete.edges
#' @param object igraph/mgraph object
#' @param es vector: edge ids (integer) or names (character)
#' @return igraph/mgraph object
#' @author ZG Zhao
#' @export
setGeneric("esdelete", function(object, es) standardGeneric("esdelete"))
#' @export
delete.edges <- function(...) esdelete(...)
#' @export
delete_edges <- function(...) esdelete(...)

setMethod("esdelete", "igraph", function(object, es) {
    igraph::delete.edges(object, es)
})

setMethod("esdelete", "xgraph", function(object, es) {
    ss1 <- 1:ecount(object) %in% es
    ss2 <- enames(object) %in% es
    ss3 <- rnames(object) %in% es
    g <- igraph::delete.edges(object, which(ss1 | ss2 | ss3))
    attributes(g) <- attributes(object)
    g
})
