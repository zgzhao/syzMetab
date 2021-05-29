#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-
# File: stat.struct.R
# Description: structural robustness/probality calculations.
# AUTHOR: ZG Zhao; zgzhao@foxmail.com
# 2020-09-01 10:40:38

XCHEM1 <- "PRI_SUBSTRATE"
XCHEM2 <- "END_PRODUCT"

.transNames <- function(chems, g) {
    xnames <- Aliases(g)
    if(! is.empty(xnames)) {
        for (nn in names(xnames)) {
            ss <- chems %in% xnames[[nn]]
            chems[ss] <- nn
        }
    }
    intersect(Chemicals(g), chems)
}

#' @export
OP_chem <- function(g, chem, substrates, ex.genes=NULL) {
    chem <- .transNames(chem, g)
    substrates <- .transNames(substrates, g)
    if(chem %in% substrates) {
        pp <- vdata(g, "PP", chem)
        return(pp)
    }
    g <- mgraph_clean(g, substrates, chem)
    from <- unlist(graph::inEdges(chem, g))
    pps <- NULL
    ## randomized order
    for(vs in sample(from)) {
        ppx <- OP_react(g, vs, chem, substrates, ex.genes)
        pps <- c(pps, ppx)
    }
    names(pps) <- NULL
    return(1 - prod(1 - pps))
}

## use with OP_chem.
## FIXME: take genes into account
OP_react <- function(g, from, to, substrates, ex.genes=NULL) {
    ename <- paste(from, to, sep="|")
    genes <- edata(g, "gene", ename)
    if(is.empty(genes))  epx <- edata(g, "EP", ename)
    else {
        genes <- setdiff(genes, ex.genes)
        ex.genes <- setdiff(union(ex.genes, genes), c(NA, "auto"))
        epx <- if(all(genes == "auto")) EPauto(g) else 1
        nn <- length(setdiff(genes, "auto"))
        epx <- epx * edgeDataDefaults(g, "EP")^nn
    }
    ppx <- OP_chem(g, from, substrates, ex.genes)
    ## 如果没有基因已排除，结果正常
    ## 如果全部基因都排除：nn=0, epx=1，传递给OP_chem的值（ppx）为0
    res <- ppx * (1 - epx)
    names(res) <- NULL
    res
}

## Convert metabolic graph (graphMET) to PGM (probability graph model).
## TODO: write a more flexible function according to source codes of gRain.
mgraph_2PGM <- function(g){
    vns <- nodes(g)
    dd <- graph::degree(g)$inDegree
    yn <- c("y", "n")
    ## vertices without parents
    pp <- list()
    for(vx in vns[dd == 0]){
        lbx <- vns[which(vns == vx)]
        pp <- c(pp, list(cptable(lbx, values = c(1, 0), levels = yn)))
    }
    
    ## causal effects
    ess <- edges(reverseEdgeDirections(g))
    for(tt in names(which(dd > 0))){
        ss <- ess[[tt]]
        if(length(ss) < 1) next
        
        dpx <- NULL
        for(i in 1:length(ss)) dpx <- c(dpx, edgesEP(g, ss[i], tt))
        lbx <- c(tt, ss)
        
        mm <- matrix(rep(c(TRUE,FALSE), length(dpx)), nrow = 2)
        mm <- expand.grid(as.data.frame(mm))
        mm <- mm[, ncol(mm):1]
        if(length(dpx) < 2) mm <- matrix(mm, nrow = 2)
        pxx <- apply(mm, 1, FUN = function(x){
            xx <- dpx[x]
            if(length(xx) < 1) rex <- 1
            else rex <- prod(xx)
            rex
        })
        pxx <- unlist(pxx)
        pxx <- as.vector(rbind(1 - pxx, pxx))
        pp <- c(pp, list(cptable(lbx, values = pxx, levels = yn)))
    }
    compileCPT(pp)    
}

#' FP (Failure Probability) is the probality of a pathway unable to get final products from given substrates. FP is calculated by probability propagation different algorithms. If more than one substrates (product) are given, they are treated equally and merged into a single "virtual substrate (product)"; all redundant nodes and edges will be removed accordingly. Please refer to \code{\link{mgraph_set_chems}} for more detail about merging of chemicals.
#'
#' When error occurs, an integer (error code) will be returned, allowing smooth running of the function in parallel computing. The codes are:
#' - 100: non graphMET object was feed or generated
#' - 200: no node found
#' - 300: no edge found
#' - 400: substrates or products not presented in the pathway
#' - 500: substrates and products not connected at all
#' @title Calculate FP of a metabolic pathway
#' @param g graphMET object
#' @param substrates character vector. Node names representing substrates or starting points of the pathway.
#' @param products character vector. Node names representing products or end points of the pathway.
#' @param simplify TRUE/FALSE. If TRUE (default), simplify network before PGM construction and FP calculation.
#' @param algorithm character-string, one of "default", "forward" and "PGM"
#' @param PPsubs numeric vector, inital PP for substrates.
#' @return
#' @author ZG Zhao
#' @export
FP_path <- function(g, substrates, products, simplify = TRUE, algorithm = "default",
                    PPsubs = 1, silent=FALSE){
    if(! silent) {
        stime <- Sys.time()
        cat("Probility is between 0 and 1. Number greater than 1 denotes error.\n")
    }

    if(! is.mgraph(g) || is.empty(g))
        return(100)

    if(is.chemSet(g)) {
        if(! (missing(substrates) || missing(products)))
           if(! silent) warning("Chemicals already set. Use old settings.")
    } else {
        g <- mgraph_set_chems(g, substrates, products)
        if(is.empty(g)) return(100)
    }

    ## Following codes handle graphMET having merged chemicals!
    g <- mgraph_clean(g)
    if(simplify) g <- mgraph_simplify(g)
    else g <- mgraph_update(g)
    if(is.empty(g)) return(100)

    ## no path from substrate to product
    if(! is.linked(g, VCHEM1, VCHEM2)) return(200)
    vdata(g, "PP") <- NA
    vdata(g, "PP", VCHEM1) <- PPsubs

    if(algorithm == "PGM") {
        ## probability graph model
        if(! gRbase::is_dag(g)) stop("The graph object is not a DAG!")
        plst <- mgraph_2PGM(g)
        chem <- names(which(graph::degree(g)$outDegree == 0))[1]
        fp <- querygrain(grain(plst), nodes = chem, type = "joint")
        fp <- fp["n"]
        names(fp) <- NULL
    } else if (algorithm == "forward"){
        ## forward P propagation
        if(! gRbase::is_dag(g)) stop("The graph object is not a DAG!")
        g <- PP_forward(g, VCHEM1)
        fp <- vdata(g, "PP", VCHEM2)
        fp <- 1 - unlist(fp)
    } else if(algorithm == "test") {
        ## test of algorithm without using simple paths saved in graph object
        ppc <- OX_chem(g, VCHEM2, VCHEM1)
        fp <- 1 - ppc
    } else {
        ppc <- OP_chem(g, VCHEM2, VCHEM1)
        fp <- 1 - ppc
    }
    names(fp) <- NULL
    if(! silent) print(Sys.time() - stime)
    fp
}

