#' plot mgraph
#'
#' Specified igraph plot function for mgraph
#' @title plot mgraph
#' @aliases plot.ggraph
#' @param g mgraph or igraph object
#' @param show.name logi. Show chemical names if TRUE (default).
#' @param ... pars passed to plot.igraph
#' @return NULL
#' @author zhao
#' @export
plot.mgraph <- function(g, show.name=TRUE,  ...) {
    options(warn=FALSE)
    opar <- par("mar")
    ipar <- igraph.options()
    par(mar=rep(0,4))
    if(show.name) {
        vcolor <- NA
        vlcol <- "black"
    } else {
        vcolor <- "#6495ED"
        vlcol <- "transparent"
    }
    if(!is.empty(Substrates(g))) {
        ss <- vnames(g) %in% Substrates(g)
        vlcol[ss] <- "red"
    }
    if(!is.empty(Products(g))) {
        ss <- vnames(g) %in% Products(g)
        vlcol[ss] <- "blue"
    }
    igraph.options(vertex.size=10,
                   vertex.color=vcolor,
                   vertex.frame.color=vcolor,
                   vertex.label.color=vlcol,
                   vertex.label.cex=1,
                   plot.layout=layout.kamada.kawai,
                   plot.margin=0,
                   edge.arrow.size = 0.5,
                   edge.arrow.width = 0.6, edge.width = 1,
                   edge.label.cex=0.8,
                   edge.label.color="red",
                   edge.color="gray40")
    plot.igraph(g, ...)
    par(mar=opar)
}

#' @export
plot.rgraph <- function(...) plot.mgraph(...)
#' @export
plot.ggraph <- function(...) plot.mgraph(...)

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

