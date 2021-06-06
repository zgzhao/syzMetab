#' Set chemicals and subset a xgraph (mgraph, rgraph or ggraph) object.
#'
#' The function do these jobs:
#' - for rgraph (reaction graph) and ggraph (gene graph), S and P are added to the graph as nodes.
#' - all edges to S or from P are deleted, and then nodes not connecting to S or P are delete.
#' - add "substrates" and "products" attributes are the graph.
#' @title set substrates and products for a graph
#' @param object xgraph object
#' @param s character vecter, names of primary substrates
#' @param p character vecter, names of end products
#' @param clean TRUE/FALSE(default). Remove futile nodes and edges if TRUE.
#' @return xgraph object
#' @author ZG Zhao
#' @export
setGeneric("xgraph_setchems", function(object, s, p, clean=FALSE) standardGeneric("xgraph_setchems"))
setMethod("xgraph_setchems", "mgraph", function(object, s, p, clean){
    if(is.chemset(object)) stop("Chemicals set already!")
    vns <- vnames(object)
    s <- intersect(s, vns)
    p <- intersect(p, vns)
    if(is.empty(s) || is.empty(p)) stop("Invalid chemical names!")
    if(! vis_connected(object, s, p))
        stop("Sources and targets are not connected in the network!")
    g <- .setChems(object, s, p, clean)
    g
})
setMethod("xgraph_setchems", "rgraph", function(object, s, p, clean){
    if(is.chemset(object)) stop("Chemicals set already!")
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    p <- intersect(p, vns)
    if(is.empty(s) || is.empty(p)) stop("Invalid chemical names!")
    if(! vis_connected(object, s, p))
        stop("Sources and targets are not connected in the network!")
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
    g <- .setChems(g, s, p, clean)
    g
})
setMethod("xgraph_setchems", "ggraph", function(object, s, p, clean){
    if(is.chemset(object)) stop("Chemicals set already!")
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    p <- intersect(p, vns)
    if(is.empty(s) || is.empty(p)) stop("Invalid chemical names!")
    if(! vis_connected(object, s, p))
        stop("Sources and targets are not connected in the network!")
    vss <- unique(c(s, p))
    g <- add.vertices(object, length(vss), name=vss)
    ## add edges according to reaction info
    rtns <- Reactions(object)
    for(rx in rtns) {
        xx <- s %in% rx[["substrate"]]
        g <- add.edges(g, expand.grid(s[xx], rx[["gene"]]))
        xx <- s %in% rx[["product"]]
        g <- add.edges(g, expand.grid(rx[["gene"]], s[xx]))
    }
    attributes(g) <- attributes(object)
    g <- .setChems(g, s, p, clean)
    g
})

.setChems <- function(g, s, p, clean) {
    Substrates(g) <- s
    Products(g) <- p
    if(! clean) return(g)

    mc.cores <- detectCores() - 1
    spp <- all_spaths_list(g, s, p, mc.cores)
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
    g
}

