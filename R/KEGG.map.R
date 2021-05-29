
#' Get KOG to gene mapping list for an organism. KEGG brite file will be download if needed.
#'
#' Refer to \code{\link{kogs_table}} if you want gene to KOG mapping table.
#' @title KOG to gene mapping list for an organism
#' @param org character, KEGG organism abbreviation
#' @param d.path file path for \code{\link{KEGG_get}}
#' @param KOGs character vector, names of KOGs. Subset if not empty.
#' @return list of genes named by KOGs.
#' @author ZG Zhao
#' @export
kogs_list <- function(org, d.path="KEGG", KOGs=NA) {
    if(is.empty(org)) return(NULL)
    if("KGDF" %in% class(org)) dtx <- org
    else  dtx <- kogs_table(org, d.path, KOGs)

    if(is.empty(dtx)) return(NULL)
    dtx <- aggregate(gene~KOG, data = dtx, FUN = unique)
    if(nrow(dtx) < 1) return(list())
    genes <- dtx$gene
    if(! is.list(genes)) genes <- as.list(genes)
    names(genes) <- dtx$KOG
    genes
}

#' Get gene to KOG mapping table for an organism. KEGG brite file will be download if needed.
#'
#' If org="ko", there are only two types of data: KOG (KEGG gene entries) and description. To provide a general interface, these data is also named "gene" and are duplicated in "ko" column. KEGG brite file for the organism will be downloaded automatically if not exists.
#' @title Gene to KOG mapping table for an organism
#' @param org character, KEGG organism abbreviation
#' @param d.path file path for \code{\link{KEGG_get}}
#' @param KOGs character vector, names of KOGs. Subset if not empty.
#' @return  data.frame with three columns: gene, KOG and desc (functional description). NOTE: Some KEGG brite file may have no gene info, NULL is return without warning!
#' @author ZG Zhao
#' @export
kogs_table <- function(org, d.path="KEGG", KOGs=NA){
    f.brite <- KEGG_get(org, d.path)
    if(! file.exists(f.brite)) stop("KEGG brite file not exists and downloaded failed!")
    kinfo <- grep("^D", readLines(f.brite), value = TRUE)
    kinfo <- gsub("^D\\s+", "", kinfo)
    if(length(kinfo) < 1) return(NULL)

    isGeneral <- grepl("^K", kinfo[1])
    results <- sapply(kinfo, FUN = function(x){
        gene <- sub("^(\\S+)\\s.+$", "\\1", x)
        if(isGeneral) {
            kog <- gene
            desc <- sub("^\\S+\\s+(.+)$", "\\1", x)
            desc <- sub("\\s+\\[.+$", "", desc)
        } else {
            kog <- NA
            desc <- NA
            if(grepl("^.+\\bK[0-9]{5}\\b.+$", x[1])) {
                kog <- sub("^.+\\b(K[0-9]{5})\\b.+$", "\\1", x)
                desc <- sub("^\\S+\\b+(.+)\\s+K[0-9]{5}.*$", "\\1", x)
            }
        }
        c(gene, kog, desc)
    })
    results <- as.data.frame(t(results))
    rownames(results) <- NULL
    colnames(results) <- c("gene", "KOG", "desc")
    results <- results[!is.na(results$KOG), ]
    KOGs <- intersect(KOGs, results$KOG)
    if(length(KOGs) > 0) results <- results[results$KOG %in% KOGs, ]
    class(results) <- c(class(results), "KGDF")
    results
}

#' Get KEGG id mapping data (gene id to uniprot or ncbi-proteinid) of a organism.
#'
#' Retrieve id mapping data from KEGG website and save locally for later use.
#' @title Function KEGG_gid2pid
#' @param org KEGG organism name (abbreviation).
#' @param d.type character: uniprot or ncbi
#' @param d.path Destination path (folder) to save or load downloaded file.
#' @return data frame with two columns: gene.id and protein.id
#' @author ZG Zhao
#' @export
KEGG_gid2pid <- function(org, d.type=c("uniprot", "ncbi"), d.path = "KEGG/maps") {
    d.type <- d.type[1]
    if(! d.type %in% c("uniprot", "ncbi"))
        stop("Allow data type of uniprot or ncbi only.")
    d.path <- file.path(d.path, d.type)
    if(d.type == "ncbi") d.type <- "ncbi-proteinid"
    
    ff <- file.path(d.path, paste0(org, ".csv"))
    if(file.exists(ff)) {
        info <- read.csv(ff)
        info[[1]] <- as.character(info[[1]])
        return(info)
    }
    if(! dir.exists(d.path)) dir.create(d.path, recursive = TRUE)

    info <- NULL
    try(info <- KEGGREST::keggConv(d.type, org, querySize = 100), silent = TRUE)
    if(is.null(info)) {
        warning(paste("Gene info retrieve failed for organism", org))
        return(invisible(NULL))
    }
    ggid <- names(info)
    names(info) <- NULL
    ppid <- sub("^.+:(.+)$", "\\1", info)
    ggid <- sub("^.+:(.+)$", "\\1", ggid)
    info <- data.frame(gene.id = ggid, protein.id = ppid)
    write.csv(info, ff, row.names = FALSE)
    info
}

#' Map uniprot ids to KEGG (K gene)
#'
#' Map uniprot ids to K genes. Databases of there species are included: ath, osa and hsa. More databases (files) can be included by db.files whose first two columns should be uniprot id and kgene.
#' @title Function uniprot2kgene
#' @param uids chr vector
#' @param db.files additional database (uniprot to k gene mapping) file.
#' @return data.frame with two columns: uniprot and kgene
#' @author ZG Zhao
#' @export
uniprot2kgene <- function(uids, db.files=NULL) {
    fdir <- file.path(path.package("zAnno"), "uniprot")
    res <- NULL
    for(ss in c("ath", "osa", "hsa")) {
        dx <- read.csv(file.path(fdir, paste0("mapped.", ss, ".csv")))
        dx <- dx[dx$uniprot %in% uids, c("uniprot", "gene")]
        dx$gene <- paste(ss, dx$gene, sep=":")
        res <- rbind(res, dx)
    }
    colnames(res) <- c("uniprot", "kgene")
    xids <- setdiff(uids, res$uniprot)
    if(length(xids) > 0 && !is.null(db.files)) {
        for (ff in db.files) {
            dx <- read.csv(ff)[, 1:2]
            colnames(dx) <- c("uniprot", "kgene")
            res <- rbind(res, dx[dx$uniprot %in% xids, ])
        }
    }
    res
}

#' Update the uniprot to kgene mapping file.
#'
#' Update the uniprot to kgene mapping file.
#' @title Function updateDB_u2k
#' @param db.file 
#' @param uids chr vector
#' @param overwirte TRUE/FALSE, write back to db.file if TRUE.
#' @return data.frame
#' @author ZG Zhao
#' @export
updateDB_u2k <- function(db.file, uids, overwirte=FALSE) {
    ## exclude proteins in these species
    fdir <- file.path(path.package("zAnno"), "uniprot")
    for(ss in c("ath", "osa", "hsa")) {
        xdb <- read.csv(file.path(fdir, paste0("mapped.", ss, ".csv")))
        uids <- setdiff(uids, xdb$uniprot)
    }

    ## exclude existing ids
    xdb <- read.csv(db.file)[, 1:2]
    colnames(xdb) <- c("uniprot", "kgene")
    uids <- setdiff(uids, xdb$uniprot)
    if(length(uids) > 0) {
        xids <- KEGGREST::keggConv("genes", paste0("uniprot:", uids), querySize = 100)
        dt <- data.frame(uniprot=names(xids), kgene=xids)
        dt$uniprot <- sub("up:", "", dt$uniprot)
        dx <- data.frame(uniprot=uids)
        dt <- merge(dt, dx, by="uniprot", all=TRUE)
        xdb <- rbind(xdb, dt)
        xdb <- xdb[!duplicated(paste(xdb$uniprot, xdb$kgene)), ]
        if(overwirte) write.csv(xdb, db.file, row.names=FALSE)
    }
    cat("Done!\n")
    xdb
}
