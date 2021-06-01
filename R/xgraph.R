#' Set chemicals and subset a xgraph (mgraph, rgraph or ggraph) object.
#'
#' The function do these jobs:
#' - for rgraph (reaction graph) and ggraph (gene graph), S and P are added to the graph as nodes.
#' - all metabolites (nodes) having no contribution in forming of P from S are deleted, and so do the reactions (edges).
#' This is done by screening nodes which are in at least one simple pathway from S to P.
#' - add "substrates" and "products" attributes are the graph.
#' @title set substrates and products for a graph
#' @param object xgraph object
#' @param s character vecter, names of primary substrates
#' @param p character vecter, names of end products
#' @param mc.cores integer, for \code{\link{all_spath_list}}, refer to \code{\link{mclapply}} for details.
#' @param n.max numeric, for \code{\link{all_spath_list}}.
#' @return xgraph object
#' @author ZG Zhao
#' @export
setGeneric("xgraph_setchems",
           function(object, s, p, mc.cores, n.max) standardGeneric("xgraph_setchems"))

setMethod("xgraph_setchems", "mgraph",
          function(object, s, p, mc.cores, n.max=50){
    if(is.chemset(object)) stop("Chemicals set already!")
    vns <- vnames(object)
    s <- intersect(s, vns)
    p <- intersect(p, vns)
    if(is.empty(s) || is.empty(p)) stop("Invalid chemical names!")
    .setChems(object, s, p, mc.cores, n.max)
})
setMethod("xgraph_setchems", "rgraph",
          function(object, s, p, mc.cores, n.max=50){
    if(is.chemset(object)) stop("Chemicals set already!")
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    p <- intersect(p, vns)
    if(is.empty(s) || is.empty(p)) stop("Invalid chemical names!")
    vss <- c(s, p)
    nn <- length(vss)
    g <- add.vertices(object, nn, name=vss)
    ## add edges according to reaction info
    rtns <- Reactions(object)
    for(vv in names(rtns)) {
        rx <- rtns[[vv]]
        ss <- s %in% rx[["substrate"]]
        if(sum(ss) > 0) g <- add.edges(g, rbind(s[ss], vv))
        ss <- p %in% rx[["product"]]
        if(sum(ss) > 0) g <- add.edges(g, rbind(vv, p[ss]))
    }
    attributes(g) <- attributes(object)
    g <- .setChems(g, s, p, mc.cores, n.max)
    g
})
setMethod("xgraph_setchems", "ggraph",
          function(object, s, p, mc.cores, n.max=50){
    if(is.chemset(object)) stop("Chemicals set already!")
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    p <- intersect(p, vns)
    if(is.empty(s) || is.empty(p)) stop("Invalid chemical names!")
    vss <- unique(c(s, p))
    g <- add.vertices(object, length(vss), name=vss)
    ## add edges according to reaction info
    rtns <- Reactions(object)
    for(rx in rtns) {
        xx <- s %in% rx[["substrate"]]
        g <- xaddEdges(g, s[xx], rx[["gene"]])
        xx <- s %in% rx[["product"]]
        g <- xaddEdges(g, rx[["gene"]], s[xx])
    }
    attributes(g) <- attributes(object)
    g <- .setChems(g, s, p, mc.cores, n.max)
    g
})


.setChems <- function(g, s, p, mc.cores, n.max){
    if(missing(mc.cores)) mc.cores <- detectCores() - 1
    if(mc.cores < 1) mc.cores <- 1

    spp <- all_spaths_list(g, s, p, mc.cores, n.max)
    vss <- all_spaths_nodes(spp)
    ess <- all_spaths_edges(spp)
    ## delete nodes
    vxx <- setdiff(vnames(g), vss)
    g <- delete.vertices(g, vxx)
    ## delete edges
    exx <- setdiff(enames(g), ess)
    g <- delete.edges(g, exx)

    ## save time-consuming results
    attr(g, "spaths") <- spp
    Substrates(g) <- s
    Products(g) <- p
    g
}
