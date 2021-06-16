#' API for downloding KEGG info with KEGG offical rest API.
#'
#' Download XML, image (png), AA sequence or NT sequence file from KEGG website. 
#' An alternative function is keggGet implemented in "KEGGREST" package.
#' @title KEGG download API
#' @param qid KEGG query id such as KEGG ortholog id like "ko00010", gene name like "ath:AT1G00010" or organism name like "ath" (for KEGG brite file only).
#' @param d.path Destination path (folder) to save or load downloaded file.
#' @param f.type Character string, file type of KEGG entry: "kgml", "aaseq", "ntseq" or "image". 
#' @param force.download TRUE/FALSE (default). Force update if TRUE and file exists.
#' @return path to the saved file.
#' @author ZG Zhao
#' @export
KEGG_get <- function(qid, d.path = "KEGG", f.type = c("kgml", "image", "htext"),
                     force.download = FALSE) {
    qid <- tolower(qid[1])
    if(grepl("^[a-z]+$", qid)) {
        ## for kegg brite file
        d.path <- .setKEGGLocalPath(d.path, "brite")
        f.path <- file.path(d.path, paste0(qid, "00001.txt"))
        xurl <- paste0("http://rest.kegg.jp/get/br:", qid, "00001")
    } else {
        ## others
        f.type <- intersect(c("kgml", "image", "htext"), tolower(f.type[1]))
        d.path <- .setKEGGLocalPath(d.path, f.type)
        if(length(f.type) < 1) return(invisible(NULL))
        
        ent <- tolower(qid)
        ext <- c(".xml", ".png", ".txt")
        names(ext) <- c("kgml", "image", "htext")
        ext <- ext[f.type]
        f.path <- file.path(d.path, paste0(ent, ext))
        
        kapi <- "http://rest.kegg.jp/get/"
        xtt <- if(f.type %in% c("kgml", "image")) paste0("/", f.type, "/") else ""
        xurl <- paste0(kapi, ent, xtt)
    }

    ## commit download
    if(! file.exists(f.path) || force.download) {
        if(! dir.exists(d.path)) dir.create(d.path, recursive = TRUE)
        try(download.file(xurl, f.path, method = "auto"), silent = TRUE, outFile = "kegg.log")
    }
    ## check after download
    if(! file.exists(f.path) || file.size(f.path) < 10) {
        warning("File does not exist or downloaded failed.")
        invisible(NULL)
    }
    return(f.path)
}


## parse KEGG XML tree infos, for internal use only.
.parseEntryList <- function(kdata) {
    itype <- xml_name(kdata)
    items <- kdata[itype == "entry"]
    if(length(items) < 1) return(NA)
    rids <- xml_attr(items, "id")
    rlist <- xml_attrs(items)
    rlist <- lapply(rlist, FUN=function(x) {
        rx <- x[setdiff(names(x), c("id", "link"))]
        rx <- gsub("\\b[a-z]+:", "", rx, ignore.case=TRUE)
        rx <- sapply(rx, FUN=function(z) sort(strsplit(z, " +")[[1]]))
        as.list(rx)
    })
    names(rlist) <- rids
    class(rlist) <- "EntryList"
    rlist
}

.parseReactionList <- function(kdata, entries) {
    itype <- xml_name(kdata)
    items <- kdata[itype == "reaction"]
    if(length(items) < 1) return(NA)
    rids <- xml_attr(items, "id")
    rlist <- lapply(items, FUN=function(x){
        xid <- xml_attr(x, "id")
        rtype <- xml_attr(x, "type") == "reversible"
        chem1 <- xml_attr(xml_find_all(x, ".//substrate"), "id")
        chem2 <- xml_attr(xml_find_all(x, ".//product"), "id")
        chem1 <- lapply(chem1, FUN=function(aa) entries[[aa]][["name"]])
        chem2 <- lapply(chem2, FUN=function(aa) entries[[aa]][["name"]])
        chem1 <- unlist(chem1)
        chem2 <- unlist(chem2)
        names(chem1) <- NULL
        names(chem2) <- NULL
        r.names <- entries[[xid]][["reaction"]]
        genes <- entries[[xid]][["name"]]
        list(substrate=chem1, product=chem2, reversible=rtype, gene=genes, name=r.names)
    })
    names(rlist) <- rids
    class(rlist) <- "ReactionList"
    attr(rlist, "raw") <- TRUE
    rlist
}

.parseGraphicsList <- function(kdata) {
    items <- xml_find_all(kdata, ".//graphics")
    ids <- xml_attr(xml_parent(items), "id")
    rlist <- xml_attrs(items)
    names(rlist) <- ids
    class(rlist) <- "GraphicList"
    rlist
}

.parseRelationList <- function(kdata) {
    itype <- xml_name(kdata)
    items <- kdata[itype == "relation"]
    if(length(items) < 1) return(NA)
    rlist <- lapply(items, FUN=function(x){
        rx1 <- xml_attrs(x)
        rx2 <- unlist(xml_attrs(xml_children(x)))
        c(rx1, rx2)
    })
    class(rlist) <- "RelationList"
    rlist
}

.filterEList <- function(entries, by, vals) {
    by <- unlist(by)[1]
    ss <- sapply(entries, FUN=function(x) {
        if(by %in% names(x)) x[by] %in% vals
        else FALSE
    })
    entries[ss]
}

## download and parse CPD info
.parseCPD <- function(KOPs, d.path="KEGG") {
    results <- NULL
    for (kpx in KOPs) {
        f.path <- KEGG_get(kpx, d.path=d.path, f.type="htext")
        dt <- readLines(f.path)
        ndx <- grep("^(COMPOUND|REFERENCE)", dt)
        ggs <- dt[ndx[1]:(ndx[2] - 1)]
        ggs <- gsub("^COMPOUND", "", ggs)
        ccs <- gsub("^ *([^ ]+).*$", "\\1", ggs)
        des <- gsub("^ *[^ ]+ +(.*)$", "\\1", ggs)
        dt <- data.frame(CPD = ccs, name = des)
        results <- rbind(results, dt)
    }
    results
}

.setKEGGLocalPath <- function(d.path, type){
    if(grepl("KEGG$", d.path, ignore.case = TRUE))
        return(file.path(d.path, type))
    ## left for user
    return(d.path)
}

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
    genes$auto <- "auto"
    genes$HLINK <- "HLINK"
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
