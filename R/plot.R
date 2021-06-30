#' plot graph
#'
#' Specified igraph plot function for mgraph
#' @title plot graph
#' @aliases plot.rgraph plot.ggraph plot.bgraph
#' @param g xgraph object
#' @param s substrates
#' @param p products
#' @param t targets, same as "p"
#' @param show.name logi. Show chemical names if TRUE (default).
#' @param gene.n TRUE/FALSE (default). show gene number on edge if TRUE.
#' @param ... pars passed to plot.igraph
#' @return NULL
#' @author zhao
#' @export
plot.mgraph <- function(g, s, t, p=t, show.name=TRUE, gene.n=FALSE, ...) {
    options(warn=FALSE)
    opar <- par("mar")
    ipar <- igraph.options()
    nv <- vcount(g)
    ne <- ecount(g)
    par(mar=rep(0,4))

    xlayout <- layout.kamada.kawai
    elty <- rep("solid", ne)
    vsize <- rep(10, nv)
    if(show.name) {
        vcolor <- rep(NA, nv)
        vlcol <- rep("black", nv)
    } else {
        vcolor <- rep("#6495ED", nv)
        vlcol <- rep("transparent", nv)
    }
    if(is.bgraph(g)) {
        xfac <- V(g)$isfac
        vsize[xfac] <- 16
        vcolor[xfac] <- "pink"
        xlayout <- layout_with_kk
    } else if(is.mgraph(g)){
        rtns <- Reactions(g)
        rtns <- rtns[rnames(g)]
        ngene <- sapply(rtns, FUN=function(x) {
            ss <- grepl("^auto", x[["gene"]])
            sum(! ss)
        })
        if(gene.n) E(g)$label <- ngene
        elty[ngene < 1] <- "dotted"
    }
    if(! missing(s)) Substrates(g) <- s
    if(! missing(p)) Products(g) <- p
    if(!is.empty(Substrates(g))) {
        ss <- vnames(g) %in% Substrates(g)
        if(show.name) vlcol[ss] <- "red"
        else vcolor[ss] <- "red"
    }
    if(!is.empty(Products(g))) {
        ss <- vnames(g) %in% Products(g)
        if(show.name) vlcol[ss] <- "blue"
        else vcolor[ss] <- "darkgreen"
    }

    igraph.options(vertex.size=vsize,
                   vertex.color=vcolor,
                   vertex.frame.color=vcolor,
                   vertex.label.color=vlcol,
                   vertex.label.cex=1,
                   plot.layout=xlayout,
                   plot.margin=0,
                   edge.arrow.size = 0.5,
                   edge.arrow.width = 1, edge.width = 2,
                   edge.label.cex=1,
                   edge.label.color="black",
                   edge.lty=elty,
                   edge.color="gray40")
    plot.igraph(g, ...)
    par(mar=opar)
}

#' @export
plot.rgraph <- function(...) plot.mgraph(...)
#' @export
plot.ggraph <- function(...) plot.mgraph(...)
#' @export
plot.bgraph <- function(...) plot.mgraph(...)


#' This function is not relevant to KEGG pathway. It is designed for illustration only.
#'
#' Install "diagram" package manually before using this function.
#' @title show reaction diagram
#' @param substrate substrate names
#' @param product product names
#' @param r.name reaction name
#' @param ... pars passed to \code{\link{diagram::plotmat}}
#' @return NULL
#' @author zhao
#' @export
plotReaction <- function(substrates, products, r.name,
                         shape.chem="ellipse", shape.react="rect",
                         lcol="gray", arr.pos=0.7, arr.col="gray", arr.tcol="transparent",
                         box.prop=0.5, shadow.size=0, relsize=0.85,
                         ...) {
    opar <- par("mar")
    par(mar=rep(0,4))
    ss <- unique(substrates)
    pp <- unique(products)
    nx <- r.name
    xnames <- c(ss, nx, pp)
    n1 <- length(ss)
    n2 <- length(pp)
    nn <- length(xnames)
    mm <- matrix(nrow=nn, ncol=nn, data=0)
    colnames(mm) <- rownames(mm) <- xnames
    mm[nx, ss] <- 1
    mm[pp, nx] <- 1
    x <- c(rep(0.1, n1), 0.5, rep(0.9, n2))
    gap <- 0.8/max(n1-1, n2-1)
    y1 <- if(n1>1) seq((1-gap*(n1-1))/2, by=gap, length=n1) else 0.5
    y2 <- if(n2>1) seq((1-gap*(n2-1))/2, by=gap, length=n2) else 0.5
    y <- c(rev(y1), 0.5, rev(y2))
    pos <- cbind(x, y)
    bt <- rep(shape.chem, nn)
    bt[xnames==nx] <- shape.react

    diagram::plotmat(mm, pos=pos, name=xnames, curve=0, box.type=bt,
                     box.prop=box.prop, shadow.size = shadow.size,
                     relsize = relsize, lcol=lcol,
                     arr.col = arr.col, arr.tcol = arr.tcol,
                     arr.pos= arr.pos, endhead=TRUE, ...)
    par(mar=opar)
}

