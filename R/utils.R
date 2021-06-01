
#' simple paths info: list, nodes and edges
#'
#' Wrapers for \code{\link{igraph::all_simple_paths}}. Parallel computing may be evoked if source or target number is more than one. In these series of functions, `all_spath_list` is the basic function; `all_spath_nodes` and `all_spath_edges` extract infos from the results of `all_spath_list` implicitly or explicitly.
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
#' @param n.max numeric, max number of nodes allowed in the graph. Don't set big number (default 50) since simple path searching for big network is very very time-consuming.
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
#' @export
all_spaths_list <- function(obj, from, to, mc.cores, n.max=50){
    if(vcount(obj) > n.max)
        stop("Node number exceeds `n.max` setting.")
    stime <- Sys.time()
    vss <- vnames(obj)
    from <- intersect(from, vss)
    to <- intersect(to, vss)
    if(length(from) < 1 || length(to) < 1) return(NULL)
    nsets <- expand.grid(from, to)
    nsets <- as.data.frame(t(nsets))
    colnames(nsets) <- NULL
    spp <- .xlapply(nsets, FUN=function(xx){
        xpp <- all_simple_paths(obj, xx[1], xx[2])
        lapply(xpp, names)
    }, mc.cores=mc.cores)
    spp <- unlist(spp, recursive = FALSE)
    xchems <- c(from, to)
    ss <- sapply(spp, FUN=function(x) sum(x %in% xchems) == 2)
    spp <- spp[ss]
    class(spp) <- c("sp.list", class(spp))
    cat("======Simpel paths search==========\n",
        "Found", length(spp), "simple paths in",
        as.numeric(Sys.time() - stime), "seconds.\n")
    spp
}

#' @export
all_spaths_edges <- function(obj, from, to, mc.cores, n.max=50) {
    if(is.xgraph(obj))
        splist <- all_spaths_list(obj, from, to, mc.cores, n.max)
    else splist <- obj

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
all_spaths_nodes <- function(obj, from, to, mc.cores, n.max=50) {
    if(is.xgraph(obj))
        splist <- all_spaths_list(obj, from, to, mc.cores, n.max)
    else splist <- obj
    vss <- unique(unlist(splist))
    names(vss) <- NULL
    class(vss) <- c("sp.nodes", class(vss))
    vss
}

.xlapply <- function(lst, FUN, mc.cores, ...) {
    if(length(lst) < 2) return(lapply(lst, FUN=FUN, ...))
    if(missing(mc.cores)) mc.cores <- detectCores() - 1
    mc.cores <- min(length(lst), mc.cores)
    mclapply(lst, FUN=FUN, mc.cores=mc.cores, ...)
}


#' This function list all gene cut sets. Remove any set of genes will disconnet primary substrates and end products.
#'
#' Set substrates and products before using this function.
#' @title list all gene cut sets
#' @param object xgraph object
#' @return list
#' - min.number
#' - min.sets
#' - all.sets
#' @author ZG Zhao
#' @export
setGeneric("all_cut_genes", function(object) standardGeneric("all_cut_genes"))
setMethod("all_cut_genes", "xgraph", function(object){
    if(! is.chemset(object)) stop("Chemicals not set yet.")
    cutlist <- sapply(min_st_separators(object), names)
    ## all simple paths from S to P should have at leat one of the cut chems
    allspp <- attr(object, "spaths")
    sels <- sapply(cutlist, FUN=function(x){
        ss <- sapply(allspp, FUN=function(y) any(x %in% y))
        all(ss)
    })
    cutlist <- cutlist[sels]

    ## map compounds to genes
    names(cutlist) <- sapply(cutlist, paste, collapse=" ")
    cutnames <- unique(unlist(cutlist))
    rlist <- Reactions(object)
    gmaps <- lapply(cutnames, FUN=function(cutchem){
        gene.in <- lapply(rlist, FUN=function(rr){
            if(any(cutchem %in% rr[["substrate"]])) return(rr[["gene"]])
            else return(NULL)
        })
        gene.in <- unique(unlist(gene.in))
        gene.out <- lapply(rlist, FUN=function(rr){
            if(any(cutchem %in% rr[["product"]])) return(rr[["gene"]])
            else return(NULL)
        })
        gene.out <- unique(unlist(gene.out))
        list(a=gene.in, b=gene.out)
    })
    names(gmaps) <- cutnames
    results <- lapply(cutlist, FUN=function(cutvector){
        nn <- length(cutvector)
        if(nn == 1) return(gmaps[[cutvector]])
        dd <- list()
        for(i in 1:nn) dd <- c(dd, list(c("a", "b")))
        dd <- as.data.frame(t(expand.grid(dd)))
        colnames(dd) <- NULL
        genes <- lapply(dd, FUN=function(x){
            genex <- NULL
            for(i in 1:nn) {
                rx <- gmaps[[cutvector[i]]]
                gns <- rx[[x[i]]]
                genex <- c(genex, gns)
            }
            sort(unique(genex))
        })
        names(genes) <- sapply(dd, paste, collapse="")
        genes
    })
    ## cleaning
    results <- unlist(results, recursive=FALSE)
    xname <- sapply(names(results), FUN=function(x){
        ss <- strsplit(x, "\\.")[[1]]
        s1 <- strsplit(ss[1], " ")[[1]]
        s2 <- strsplit(ss[2], "")[[1]]
        paste0(s1, s2, collapse=" ")
    })
    names(results) <- xname
    glen <- sapply(results, length)
    minn <- min(glen)
    list(min.number=minn, min.sets=results[glen == minn], all.sets=results)
})
