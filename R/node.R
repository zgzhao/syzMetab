
#' get node names
#'
#' Similar functions: \code{\link{vnames}}, \code{\link{enames}}, \code{\link{rnames}}, \code{\link{vcount}}, \code{\link{ecount}}
#' @title node names
#' @param object mgraph object
#' @author ZG Zhao
#' @export
setGeneric("vnames", function(object) standardGeneric("vnames"))
setMethod("vnames", "igraph", function(object){
    as_ids(V(object))
})
setMethod("vnames", "mgraph", function(object){
    as_ids(V(object))
})

#' Get or set node data/attributes
#'
#' Friendly function for node data/attributes retrieving or setting. See "Example".
#' @title node data
#' @param g igraph/mgraph object
#' @param a.name character, attribute name
#' @param v.names vector of character (node names) or integer (indices)
#' @seealso \code{\link{edata}}
#' @author ZG Zhao
#' @examples
#' library(gmetab)
#' d.path <- file.path(path.package("gmetab"), "KEGG")
#' gg <- mgraph_from_kos("ko00010", d.path)
#' ## igraph style
#' V(gg)$name
#' V(gg)$EP <- sample(vcount(gg))
#' V(gg)$EP[1:3]
#' ## gmetab style
#' vdata(gg, "name")
#' vdata(gg, "name", 1:3)
#' (vv <- vnames(gg))
#' vdata(gg, "EP", "C00022")
#' vdata(gg, "EP") <- 1
#' vdata(gg, "EP")
#' @export
vdata <- function(g, a.name, v.names) {
    if(! is.xgraph(g)) return(NULL)
    a.name <- a.name[1]
    ee <- substitute(V(g)$`x`, list(x=a.name))
    rx <- eval(ee)
    if(is.null(rx) || missing(v.names)) return(rx)
    ss <- (1:vcount(g) %in% v.names) | (vnames(g) %in% v.names)
    rx[v.names]
}

#' @export
`vdata<-` <- function(g, a.name, v.names, value) {
    a.name <- a.name[1]
    if(missing(v.names)) {
        ee <- substitute(V(g)$`attx` <- value, list(attx=a.name, value=value))
    } else {
        ss <- (1:vcount(g) %in% v.names) | (vnames(g) %in% v.names)
        ndx <- which(ss)
        ee <- substitute(V(g)$`attx`[n] <- value, list(attx=a.name, n=ndx, value=value))
    }
    eval(ee)
    g
}

#' delete vertices/nodes from igraph or mgraph object
#'
#' Refer to \code{\link{igraph::delete.vertices}}
#' @title delete vertices
#' @aliases delete_vertices delete.vertices
#' @param object igraph/mgraph object
#' @param vs vector: vertex ids (integer) or names (character)
#' @return igraph/mgraph object
#' @author ZG Zhao
#' @export
setGeneric("vsdelete", function(object, vs) standardGeneric("vsdelete"))
#' @export
delete.vertices <- function(...) vsdelete(...)
#' @export
delete_vertices <- function(...) vsdelete(...)
#' @export
delete.chemicals <- function(...) vsdelete(...)
#' @export
delete_chemicals <- function(...) vsdelete(...)

setMethod("vsdelete", "igraph", function(object, vs) {
    igraph::delete.vertices(object, vs)
})

setMethod("vsdelete", "mgraph", function(object, vs) {
    g <- igraph::delete.vertices(object, vs)
    attributes(g) <- attributes(object)
    g
})

#' @export
nodes_linked <- function(g, from, to) {
    nss <- nodes(g)
    from <- intersect(nss, from)
    to <- intersect(nss, to)
    if(length(from) < 1 || length(to) < 1) return(FALSE)

    g <- trans_graph(g, mode[1])
    res <- sapply(from, FUN = function(ss){
        tgs <- nodesAcc(g, ss)[[1]]
        to %in% tgs
    })
    if(! is.matrix(res)) {
        res <- matrix(res, ncol=length(from))
    }
    colnames(res) <- from
    rownames(res) <- to
    res <- t(res)
    res
}

#' @export
nodes_accessible <- function(g, from, to) {
    is.linked(g, from, to, mode)
}

#' @export
nodes_adjacent <- function(g, from, to) {
    nss <- nodes(g)
    from <- intersect(nss, from)
    to <- intersect(nss, to)
    if(length(from) < 1 || length(to) < 1) return(FALSE)

    elist <- graph::adj(trans_graph(g, mode[1]), from)
    res <- lapply(elist, FUN = function(ff) { to %in% ff })
    res <- t(as.matrix(as.data.frame(res)))
    rownames(res) <- from
    colnames(res) <- to
    res
}

#' @export
nodesAcc <- function(g, ns){
    ns <- intersect(ns, nodes(g))
    if(length(ns) < 1) return(NULL)

    rex <- graph::acc(g, ns)
    rex <- lapply(rex, names)
    rex
}

#' @export
nodes_interLinked <- function(g, from, to){
    nss <- nodes(g)
    from <- intersect(nss, from)
    to <- intersect(nss, to)
    if(length(from) < 1 || length(to) < 1) return(NULL)

    ## nodes reached by sources
    nds1 <- nodesAcc(g, from)
    ss <- sapply(nds1, function(x) any(x %in% to))
    nds1 <- nds1[ss]
    from <- names(nds1)
    to <- intersect(to, unlist(nds1))
    if(length(from) < 1 || length(to) < 1) return(NULL)

    ## nodes reach to targets
    nds2 <- nodesAcc(reverseEdgeDirections(g), to)
    nds <- intersect(unlist(nds1), unlist(nds2))
    unique(c(from, nds, to))
}

#' @export
nodesAdjacent <- function(g, vs, mode = c("out", "out", "both")) {
    elist <- graph::adj(trans_graph(g, mode[1]), vs)
    vadj <- unique(unlist(elist))
    setdiff(unique(sort(vadj)), vs)
}

#' @export
nodesChem <- function(g){
    rr <- nodeData(g, nodes(g), "ischem")
    rr <- unlist(rr)
    names(rr)[rr]
}

#' @export
nodesReaction <- function(g){
    rr <- nodeData(g, nodes(g), "ischem")
    rr <- unlist(rr)
    names(rr)[! rr]
}

.nodeGenes <- function(g, n){
    if(length(n) != 1) stop("One node name only!")
    nodeData(g, n, "genes")[[1]]
}
#' @export
nodesHub <- function(g, from, to, mode = c("out", "in", "both")) {
    vss <- nodes(g)
    from <- intersect(from, vss)
    to <- intersect(to, vss)
    if(length(from) < 1 || length(to) < 1) return(NULL)
    g <- as(g, "igraph")
    ee <- NULL
    for(ss in from) {
        spp <- igraph::all_shortest_paths(g, ss, to, mode = mode[1])$res
        if(length(spp) < 1) next
        ex <- lapply(spp, names)
        ee <- c(ee, ex)
    }
    if(length(ee) < 2) return(unlist(ee))
    res <- ee[[1]]
    for(i in 2:length(ee)) res <- intersect(res, ee[[i]])
    res <- setdiff(res, c(from, to))
    res
}

#' @export
nodesHubNearest <- function(g, v, hubs){
    if(v %in% hubs) stop("The node is already a hubs!")
    if(is.null(hubs)) stop("No hubs.")

    nn <- numNodes(g)
    res <- NULL
    vs <- v
    while(length(res) < 2) {
        vx <- nodesAdjacent(g, v, "both")
        ## no new nodes
        if(length(vs) == length(union(vs, vx))) break
        vs <- union(vs, vx)
        ## prevent search in same direction of hub got
        res <- intersect(vs, hubs)
        v <- setdiff(vs, hubs)
    }
    hubs[hubs %in% res]
}

