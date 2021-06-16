#' Subset graph with given sources (primary substrates) and targets (products)
#'
#' - In "setorg_graph", reactions without mapping genes are labeled with "absent". This function removes the "absent" links in the graph.
#' - After setting sources and targets, nodes to sources or from targets should be removed.
#' - Nodes not accessed by sources or targets should be removed also.
#' - If clean=TRUE, nodes and edges not in paths from sources to targets are further removed by calling \code{\link{all_spaths_list}}.
#' @title subset xgraph
#' @param object graph object
#' @param s source nodes
#' @param t target nodes
#' @return graph object
#' @author ZG Zhao
#' @export
#' @examples
#' ## NOT RUN
#' library(gmetab)
#' gm <- make_mgraph("ko00010")
#' chem1 <- c("C00031", "C00221", "C00267", "C01172", "C01451", "C06186")
#' chem2 <- "C00022"
#' gx <- subset_graph(gm, s=chem1, t=chem2)
#' par(mfrow=c(2,1))
#' plot(gm)
#' plot(gx)
setGeneric("subset_graph", function(object, s, t, clean=FALSE) standardGeneric("subset_graph"))
setMethod("subset_graph", "NULL", function(object, s, t, clean) return(NULL))
setMethod("subset_graph", "mgraph", function(object, s, t, clean){
    vns <- vnames(object)
    s <- intersect(s, vns)
    t <- intersect(t, vns)
    if(is.empty(s) || is.empty(t)) return(emptyGraph(object))
    if(! vis_connected(object, s, t)) return(emptyGraph(object))
    g <- cleanPrimary(object, s, t)
    if(clean) g <- cleanThorough(object, s, t)
    g
})
## TODO: if s and t are not chemicals
setMethod("subset_graph", "rgraph", function(object, s, t, clean){
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    t <- intersect(t, vns)
    if(is.empty(s) || is.empty(t)) return(emptyGraph(object))
    if(! vis_connected(object, s, t)) return(emptyGraph(object))
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
    g <- cleanPrimary(g, s, t)
    if(clean) g <- cleanThorough(object, s, t)
    g
})
## TODO
setMethod("subset_graph", "ggraph", function(object, s, t, clean){
    ## add nodes to reaction graph
    vns <- Compounds(object)
    s <- intersect(s, vns)
    t <- intersect(t, vns)
    if(is.empty(s) || is.empty(t)) return(emptyGraph(object))
    if(! vis_connected(object, s, t)) return(emptyGraph(object))
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
    g <- cleanPrimary(g, s, t)
    if(clean) g <- cleanThorough(object, s, t)
    g
})

cleanThorough <- function(object, s, t) {
    if(! vis_connected(object, s, t)) return(emptyGraph(object))
    mc.cores <- detectCores() - 1
    spp <- all_spaths_list(object, s, t)
    vss <- all_spaths_nodes(spp)
    ess <- all_spaths_edges(spp)
    ## delete nodes
    vxx <- setdiff(vnames(object), vss)
    object <- delete.vertices(object, vxx)
    ## delete edges
    exx <- setdiff(enames(object), ess)
    object <- delete.edges(object, exx)
    ## save time-consuming results
    attr(object, "spaths") <- spp
    object
}

cleanPrimary <- function(object, s, t){
    if(! vis_connected(object, s, t)) return(emptyGraph(object))
    ## remove absent reactions
    rtns <- Reactions(object)
    if(! is.empty(rtns)) {
        ss <- sapply(rtns, FUN=function(x) {
            all(x[["gene"]] %in% "absent")
        })
        rtns <- rtns[! ss]
        rnx <- names(rtns)
        rset <- Reactions(object, list.only=FALSE)
        rset@reaction <- rtns
        attr(object, "reactions") <- rset
        e.names <- enames(object)
        r.names <- rnames(object)
        exx <- e.names[! r.names %in% rnx]
        object <- delete.edges(object, exx)
    }
    ## remove not connected nodes
    vss <- vs_accessed_by(object, s, "out")
    vss <- c(vss, vs_accessed_by(object, t, "in"))
    vxx <- setdiff(vnames(object), vss)
    if(! is.empty(vxx)) object <- delete.vertices(object, vxx)
    ## remove edges to sources or from targets
    s <- intersect(s, vnames(object))
    t <- intersect(t, vnames(object))
    vxx <- vs_adjacent(object, s, "in")
    for(aa in s) {
        vtt <- vxx[[aa]]
        if(is.empty(vtt)) next
        exx <- paste(vtt, aa, sep="|")
        object <- delete.edges(object, exx)
    }
    vxx <- vs_adjacent(object, t, "out")
    for(aa in t) {
        vtt <- vxx[[aa]]
        if(is.empty(vtt)) next
        exx <- paste(aa, vtt, sep="|")
        object <- delete.edges(object, exx)
    }
    Substrates(object) <- s
    Products(object) <- t
    object
}

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
#' @examples
#' ## NOT RUN
#' library(gmetab)
#' gm <- make_mgraph("ko00010")
#' gx <- subset_graph(gm, org="ath")
#' par(mfrow=c(2,1))
#' plot(gm)
#' ## edges (reactions) with zero gene are show in dotted lines.
#' plot(gx)
setGeneric("setorg_graph", function(object, org, d.path="KEGG") standardGeneric("setorg_graph"))
#' @export
setMethod("setorg_graph", c("mgraph", "character"), function(object, org, d.path){
    if(Organism(object) != "ko") stop("Not a generic metabolic graph!")
    rtns <- Reactions(object)
    org <- org[1]
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    rtns <- lapply(rtns, FUN=function(rr){
        genes <- rr$gene
        genex <- grep("^auto.+", genes, value=TRUE)
        genes <- intersect(genes, kogs)
        if(!is.empty(genes)) genes <- unlist(gmap[genes])
        genes <- c(genes, genex)
        if(is.empty(genes)) genes <- "absent"
        names(genes) <- NULL
        rr$gene <- .setGeneNames(genes)
        rr$reversible <- FALSE
        rr
    })
    ## ensure returning result is a graph!
    ## ss <- sapply(rtns, FUN=function(x) ! is.empty(x$gene))
    ## if(sum(ss) < 1) return(emptyGraph(object))
    ## rtns <- rtns[ss]

    robj <- make_rset(rtns, org)
    gx <- make_mgraph(robj)
    if(is.chemset(object)) {
        s <- intersect(Substrates(object), vnames(gx))
        p <- intersect(Products(object), vnames(gx))
        if(is.empty(s) || is.empty(object))
            gx <- emptyGraph(gx)
        else
            gx <- subset_graph(gx, s, p)
    }
    gx
})
