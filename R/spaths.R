#' @name all_spaths_list
#' @title info of all simple paths
#' @description simple paths info: list, nodes and edges
#' @details Wrapers for \code{\link{igraph::all_simple_paths}}. Parallel computing may be evoked if source or target number is more than one. In these series of functions, `all_spaths_list` is the basic function; `all_spaths_nodes` and `all_spaths_edges` extract infos from the results of `all_spaths_list` implicitly or explicitly.
#' @aliases all_spaths_nodes all_spaths_edges
#' @usage
#' all_spaths_list(obj, from, to, mc.cores)
#' all_spaths_nodes(obj, from, to, mc.cores)
#' all_spaths_nodes(obj)
#' all_spaths_edges(obj, from, to, mc.cores)
#' all_spaths_edges(obj)
#' @param object xgraph object or list returned by all_spaths_list
#' @param from character vector of length 1, name of start node
#' @param to character vector of length 1, name of end node
#' @param mc.cores refer to \code{\link{mclapply}}
#' @param n.max numeric, max number of nodes allowed in the graph. Don't set big number (default 1000) since simple path searching for big network is very very time-consuming.
#' @return character list or vector
#' @author ZG Zhao
#' @examples
#' library(gmetab)
#' d.path <- file.path(path.package("gmetab"), "KEGG")
#' gg <- make_mgraph("ko00010", d.path)
#' vv <- vnames(gg)
#' (splist <- all_spaths_list(gg, vv[1], vv[6]))
#' all_spaths_nodes(splist)
#' all_spaths_nodes(gg, vv[1], vv[6])
#' all_spaths_edges(splist)
#' all_spaths_edges(gg, vv[1], vv[6])
NULL


#' @export
setGeneric("all_spaths_list",
           function(object, from, to, mc.cores, n.max=1000) standardGeneric("all_spaths_list"))
setMethod("all_spaths_list", "NULL", function(object, from, to, mc.cores, n.max) return(NULL))
setMethod("all_spaths_list", "xgraph",
          function(object, from, to, mc.cores, n.max){
    ## NOTE for maintaining: "clean" must be FALSE!
    object <- subset_graph(object, from, to, clean=FALSE)
    if(is.empty(object)) return(NULL)

    if(vcount(object) > n.max)
        stop("Node number exceeds `n.max` setting.")
    stime <- Sys.time()
    vss <- vnames(object)
    from <- intersect(from, vss)
    to <- intersect(to, vss)
    if(length(from) < 1 || length(to) < 1) return(NULL)
    nsets <- expand.grid(from, to)
    nsets <- as.data.frame(t(nsets))
    colnames(nsets) <- NULL
    spp <- .xlapply(nsets, FUN=function(xx){
        xpp <- all_simple_paths(object, xx[1], xx[2])
        lapply(xpp, names)
    }, mc.cores=mc.cores)
    spp <- unlist(spp, recursive = FALSE)
    xchems <- c(from, to)
    ss <- sapply(spp, FUN=function(x) sum(x %in% xchems) == 2)
    spp <- spp[ss]
    class(spp) <- c("sp.list", class(spp))
    ## cat("======Simpel paths search==========\n",
    ##     "Found", length(spp), "simple paths in",
    ##     as.numeric(Sys.time() - stime), "seconds.\n")
    spp
})

#' @export
all_spaths_edges <- function(object, from, to, mc.cores, n.max=1000) {
    if(is.xgraph(object))
        splist <- all_spaths_list(object, from, to, mc.cores, n.max)
    else splist <- object

    if(is.empty(splist)) return(NULL)
    ess <- sapply(splist, FUN=function(pp){
        nn <- length(pp)
        if(nn < 3) return(paste(pp, collapse="|"))
        paste(pp[1:(nn - 1)], pp[2:nn], sep="|")
    })
    ess <- unique(unlist(ess))
    names(ess) <- NULL
    class(ess) <- c("sp.edges", class(ess))
    ess
}

#' @export
all_spaths_nodes <- function(object, from, to, mc.cores, n.max=1000) {
    if(is.xgraph(object))
        splist <- all_spaths_list(object, from, to, mc.cores, n.max)
    else splist <- object
    vss <- unique(unlist(splist))
    names(vss) <- NULL
    class(vss) <- c("sp.nodes", class(vss))
    vss
}

