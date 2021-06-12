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
NULL


#' @export
setGeneric("all_spaths_list",
           function(object, from, to, mc.cores, n.max=50) standardGeneric("all_spaths_list"))
setMethod("all_spaths_list", "xgraph",
          function(object, from, to, mc.cores, n.max){
    object <- .preClean(object, from, to)
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
all_spaths_edges <- function(object, from, to, mc.cores, n.max=50) {
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
all_spaths_nodes <- function(object, from, to, mc.cores, n.max=50) {
    if(is.xgraph(object))
        splist <- all_spaths_list(object, from, to, mc.cores, n.max)
    else splist <- object
    vss <- unique(unlist(splist))
    names(vss) <- NULL
    class(vss) <- c("sp.nodes", class(vss))
    vss
}

#' Get st_cuts info: edges and genes
#'
#' This is a wrapper function of `st_cuts` define in igraph package.
#' @title list all gene cut sets
#' @param object xgraph object
#' @param s source node
#' @param t target node
#' @return list
#' - min.genes
#' - min.edges
#' - edge.cuts
#' - gene.cuts
#' @author ZG Zhao
#' @export
setGeneric("all_stcuts", function(object, s, t) standardGeneric("all_stcuts"))
setMethod("all_stcuts", "mgraph", function(object, s, t){
    s <- intersect(vnames(object), s)
    p <- intersect(vnames(object), t)
    if(is.empty(s) || is.empty(p)) return(NULL)
    if(! vis_connected(object, s, p)) return(NULL)

    ## TIPS: st_cuts allows one s and one p only.
    ## Merge substrates and products, respectively.
    nns <- length(s)
    nnp <- length(p)
    if(nns > 0) {
        object <- add.edges(object, expand.grid("VCHEM1", s))
        s <- "VCHEM1"
    }
    if(nnp > 0) {
        object <- add.edges(object, expand.grid(p, "VCHEM2"))
        p <- "VCHEM2"
    }
    ## finding st-cuts
    xcuts <- st_cuts(object, s, p)
    ecuts <- lapply(xcuts[[1]], as_ids)
    vcuts <- lapply(xcuts[[2]], as_ids)

    ## map edge cuts to gene cuts
    r.names <- rnames(object)
    e.names <- enames(object)
    rlist <- Reactions(object)
    gcuts <- lapply(ecuts, FUN=function(enn){
        if(grepl("VCHEM", paste(enn, collapse=" "))) return(NULL)
        ss <- r.names[e.names %in% enn]
        gns <- lapply(rlist[ss], FUN=function(rr) rr[["gene"]])
        gns <- unique(unlist(gns))
        ## cut set with autonomous reaction CANNOT be interrupted!
        if(any(c("auto", "HLINK") %in% gns)) return(NULL)
        names(gns) <- NULL
        gns
    })
    ss <- sapply(gcuts, FUN=function(x) !is.empty(x))
    if(sum(ss) < 1) return(NULL)
    ecuts <- ecuts[ss]
    gcuts <- gcuts[ss]
    names(gcuts) <- NULL
    new("stcuts", edges=ecuts, genes=gcuts)
})

## pre-cleaning network especially before simple path searching
.preClean <- function(g, s, p) {
    ## remove not connected nodes
    vss <- vs_accessed_by(g, s, "out")
    vss <- c(vss, vs_accessed_by(g, p, "in"))
    vxx <- setdiff(vnames(g), vss)
    if(! is.empty(vxx)) g <- delete.vertices(g, vxx)
    ## remove edges to sources or from targets
    s <- intersect(s, vnames(g))
    p <- intersect(p, vnames(g))
    vxx <- vs_adjacent(g, s, "in")
    for(aa in s) {
        vtt <- vxx[[aa]]
        if(is.empty(vtt)) next
        exx <- paste(vtt, aa, sep="|")
        g <- delete.edges(g, exx)
    }
    vxx <- vs_adjacent(g, p, "out")
    for(aa in p) {
        vtt <- vxx[[aa]]
        if(is.empty(vtt)) next
        exx <- paste(aa, vtt, sep="|")
        g <- delete.edges(g, exx)
    }
    g
}

