% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a_classes.R
\name{classes}
\alias{classes}
\title{package class designs}
\description{
Virtual classes designed in mgraph package.
}
\details{
Virtual classes:
\itemize{
\item EntryList:: list
\item ReactionList:: list
\item RelationList:: list
\item GraphicList:: list
\item KLists:: EntryList, ReactionList, RelationList or GraphicList
\item igraph: igraph
\item mgraph: metabolite graph
\item bgraph: bipartite graph of s-t cuts of genes
\item rgraph: reaction graph
\item ggraph: gene graph
\item xgraph: virtual class of mgraph, rgraph, ggraph and bgraph
}

S4 Classes:
\itemize{
\item KDataSet: complete set of pathway data parsed from a KEGG xml file
\item ReactionSet: reactions container used in xgraph object
}

Short names or aliases:
\itemize{
\item kdset: for \code{KDataSet}
\item klist: for \code{KLists}
\item rset: for \code{ReactionSet}
\item rdata/rnames/rcount: reaction data/names/count
\item edata/enames/ecount: edge data/names/count
\item vdata/vnames/vcount: vertex data/names/count
}
}
\examples{
library(syzMetab)
showClass("KDataSet")
showClass("ReactionSet")
showClass("mgraph")
showClass("xgraph")
}
\author{
zhao
}
