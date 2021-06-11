#' Subset graph with given sources (primary substrates) and targets (products)
#'
#' The function do these jobs:
#' - for rgraph (reaction graph) and ggraph (gene graph), S and P are added to the graph as nodes.
#' - all edges to S or from P are deleted, and then nodes not connecting to S or P are delete.
#' - add "substrates" and "products" attributes are the graph.
#' @title subset xgraph
#' @param object xgraph object
#' @param s character vecter, source ndoes
#' @param t character vecter, target nodes
#' @return xgraph object
#' @author ZG Zhao
#' @export
setGeneric("xgraph_subset", function(object, s, t) standardGeneric("xgraph_subset"))
setMethod("xgraph_subset", "mgraph", function(object, s, t){
    if(is.chemset(object)) stop("Chemicals set already!")
    vns <- vnames(object)
    s <- intersect(s, vns)
    t <- intersect(t, vns)
    if(is.empty(s) || is.empty(t)) stop("Invalid chemical names!")
    if(! vis_connected(object, s, t))
        stop("Sources and targets are not connected in the network!")
    g <- .setChems(object, s, t)
    g
})
setMethod("xgraph_subset", "rgraph", function(object, s, t){
    if(is.chemset(object)) stop("Chemicals set already!")
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    t <- intersect(t, vns)
    if(is.empty(s) || is.empty(t)) stop("Invalid chemical names!")
    if(! vis_connected(object, s, t))
        stop("Sources and targets are not connected in the network!")
    vss <- c(s, t)
    nn <- length(vss)
    g <- add.vertices(object, nn, name=vss)
    ## add edges according to reaction info
    rtns <- Reactions(object)
    for(vv in names(rtns)) {
        rx <- rtns[[vv]]
        ss <- s %in% rx[["substrate"]]
        if(sum(ss) > 0) g <- add.edges(g, rbind(s[ss], vv))
        ss <- t %in% rx[["product"]]
        if(sum(ss) > 0) g <- add.edges(g, rbind(vv, t[ss]))
    }
    attributes(g) <- attributes(object)
    g <- .setChems(g, s, t)
    g
})
setMethod("xgraph_subset", "ggraph", function(object, s, t){
    if(is.chemset(object)) stop("Chemicals set already!")
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    t <- intersect(t, vns)
    if(is.empty(s) || is.empty(t)) stop("Invalid chemical names!")
    if(! vis_connected(object, s, t))
        stop("Sources and targets are not connected in the network!")
    vss <- unique(c(s, t))
    g <- add.vertices(object, length(vss), name=vss)
    ## add edges according to reaction info
    rtns <- Reactions(object)
    for(rx in rtns) {
        xx <- s %in% rx[["substrate"]]
        g <- add.edges(g, expand.grid(s[xx], rx[["gene"]]))
        xx <- t %in% rx[["product"]]
        g <- add.edges(g, expand.grid(rx[["gene"]], t[xx]))
    }
    attributes(g) <- attributes(object)
    g <- .setChems(g, s, t)
    g
})

#' Normalize KEGG generic pathway to species specific pathway.
#'
#' Adjust and filter reactions, and rebuild graph object from reactions list.
#' @title KO-graph to species-specific graph
#' @param object xgraph object
#' @param org character
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return graph with organism set
#' @author ZG Zhao
#' @export
setGeneric("xgraph_orgset", function(object, org, d.path="KEGG") standardGeneric("xgraph_orgset"))
#' @export
mgraph_x_org <- function(...) xgraph_subset(...)
setMethod("xgraph_orgset", c("mgraph", "character"), function(object, org, d.path){
    if(Organism(object) != "ko") stop("Not a generic metabolic graph!")
    rtns <- Reactions(object)
    org <- org[1]
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    rtns <- lapply(rtns, FUN=function(rr){
        genes <- intersect(rr$gene, kogs)
        genes <- if(length(genes) < 1) NA else unlist(gmap[genes])
        names(genes) <- NULL
        rr$gene <- genes
        rr$reversible <- FALSE
        rr
    })
    ss <- sapply(rtns, FUN=function(x) ! is.empty(x$gene))
    ## ensure returning result is a graph!
    if(sum(ss) < 1) return(.emptyGraph(object))

    rtns <- rtns[ss]
    robj <- make_rset(rtns, org)
    gx <- make_mgraph(robj)
    if(is.chemset(object)) {
        s <- intersect(Substrates(object), vnames(gx))
        p <- intersect(Products(object), vnames(gx))
        if(is.empty(s) || is.empty(object))
            gx <- .emptyGraph(gx)
        else
            gx <- xgraph_subset(gx, s, p)
    }
    gx
})

.setChems <- function(g, s, t, clean) {
    Substrates(g) <- s
    Products(g) <- t
    mc.cores <- detectCores() - 1
    spp <- all_spaths_list(g, s, t, mc.cores)
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

