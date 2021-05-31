#' simple paths info: list, nodes and edges
#'
#' Friendly wrapers for \code{\link{igraph::all_simple_paths}}.
#' @title info of all simple paths
#' @aliases all_spaths_nodes all_spaths_edges
#' @usage
#' all_spath_list(obj, from, to, mc.cores)
#' all_spath_nodes(obj, from, to, mc.cores)
#' all_spath_nodes(obj)
#' all_spath_edges(obj, from, to, mc.cores)
#' all_spath_edges(obj)
#' @param obj igraph/mgraph object or list returned by all_spaths_list
#' @param from character vector of length 1, name of start node
#' @param to character vector of length 1, name of end node
#' @param mc.cores refer to \code{\link{mclapply}}
#' @return character list or vector
#' @author ZG Zhao
#' @examples
#' library(gmetab)
#' d.path <- file.path(path.package("gmetab"), "KEGG")
#' gg <- mgraph_from_kos("ko00010", d.path)
#' vv <- vnames(gg)
#' (splist <- all_spaths_list(gg, vv[1], vv[6]))
#' all_spaths_nodes(splist)
#' all_spaths_nodes(gg, vv[1], vv[6])
#' all_spaths_edges(splist)
#' all_spaths_edges(gg, vv[1], vv[6])
#' @export
all_spaths_list <- function(obj, from, to, mc.cores=detectCores() - 1){
    vss <- Chemicals(obj)
    from <- intersect(from, vss)
    to <- intersect(to, vss)
    if(length(from) < 1 || length(to) < 1) return(NULL)
    mc.cores <- min(length(from), mc.cores)
    spp <- mclapply(from, FUN=function(ss){
        xpp <- all_simple_paths(obj, ss, to)
        lapply(xpp, names)
    }, mc.cores=mc.cores)
    spp <- unlist(spp, recursive = FALSE)
    xchems <- c(from, to)
    ss <- sapply(spp, FUN=function(x) sum(x %in% xchems) == 2)
    spp[ss]
}

#' @export
all_spaths_edges <- function(obj, from, to, mc.cores=detectCores() - 1) {
    if(is.xgraph(obj))
        splist <- all_spaths_list(obj, from, to, mc.cores)
    else splist <- obj

    if(is.empty(splist)) return(NULL)
    ess <- sapply(splist, FUN=function(pp){
        nn <- length(pp)
        if(nn < 3) return(paste(pp, collapse="|"))
        paste(pp[1:(nn - 1)], pp[2:nn], sep="|")
    })
    ess <- unlist(ess)
    names(ess) <- NULL
    ess
}

#' @export
all_spaths_nodes <- function(obj, from, to, mc.cores=detectCores() - 1) {
    if(is.xgraph(obj))
        splist <- all_spaths_list(obj, from, to, mc.cores)
    else splist <- obj
    vss <- unlist(splist)
    names(vss) <- NULL
    vss
}
