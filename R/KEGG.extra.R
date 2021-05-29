## parse ko00001.keg file to data.frame
#' @export
parseKeg <- function(kegfile) {
    ## parse text
    klines <- grep("^[A-D]", readLines(kegfile), value = TRUE)
    klines <- sub("^([A-D])", "\\1 ", klines)
    klines <- grep("^[A-D] +\\w+ +\\w+", klines, value=TRUE)
    results <- t(sapply(klines, FUN = function(x){
        klevel <- sub("^([A-D]).+$", "\\1", x)
        kindex <- sub("^[A-D] +(\\w+) .*$", "\\1", x)
        kdesc <- sub("^[A-D] +\\w+ +(.*)$", "\\1", x)
        c(klevel, kindex, kdesc)
    }))
    results <- as.data.frame(results)
    colnames(results) <- c("level", "ko", "description")
    rownames(results) <- NULL
    results$description <- sub(" \\[.+$", "", results$description)
    results$ko <- toupper(results$ko)

    ## ABC hierarchical
    lvs <- results$level
    nn <- nrow(results)
    for(ll in LETTERS[1:3]){
        results[[ll]] <- NA
        ndx <- which(lvs == ll)
        for(i in 1:length(ndx)){
            ff <- ndx[i] + 1
            results[[ll]][ff:nn] <- results$ko[ndx[i]]
        }
    }
    for(i in 1:3) {
        ss <- lvs == LETTERS[i]
        for(j in i:3) results[ss, LETTERS[j]] <- NA
    }

    ## exclude unknowns
    kos <- c("09100", "09120", "09130", "09140")
    ss <- results$ko %in% kos | results$A %in% kos
    keg_na.omit(results[ss, ])
}

## optional: omit items without kgenes
#' @export
keg_na.omit <- function(kdf) {
    kos <- unlist(kdf[, c("A", "B", "C")])
    kdf[kdf$ko %in% kos | kdf$level == "D", ]
}

## convert to hierarchical list
#' @export
kdf2list <- function(df_kegg){
    kcc <- df_kegg$ko[df_kegg$level == "C"]
    kcx <- lapply(kcc, FUN=function(x){
        ss <- df_kegg$C %in% x
        df_kegg$ko[ss]
    })
    names(kcx) <- kcc
    kbb <- df_kegg$ko[df_kegg$level == "B"]
    kbx <- lapply(kbb, FUN=function(x){
        ss <- df_kegg$B %in% x
        kox <- df_kegg$ko[ss]
        kox <- intersect(kcc, kox) # C level only
        kcx[kox]
    })
    names(kbx) <- kbb
    kaa <- df_kegg$ko[df_kegg$level == "A"]
    kax <- lapply(kaa, FUN=function(x){
        ss <- df_kegg$A %in% x
        kox <- df_kegg$ko[ss]
        kox <- intersect(kbb, kox) # B level only
        kbx[kox]
    })
    names(kax) <- kaa
    kax
}

#' @export
kdf2geneList<- function(df_kegg){
    kos <- grep("^[0-9]", df_kegg$ko, value=TRUE)
    kdx <- df_kegg[df_kegg$level == "D", ]
    rex <- mclapply(kos, FUN=function(kx){
        ss <- apply(kdx, 1, FUN=function(rr) any(rr %in% kx))
        kdx$ko[ss]
    })
    names(rex) <- kos
    rex
}


#' Show or extract KEGG pathway info.
#'
#' A data.frame object "keggs" is bound with the package. koInfo is used to filter this object.
#' @title Function koInfo
#' @param kos KO indices to filer the ko column of keggs object.
#' @param keyword keyword (regular expression) to filter the description column of keggs object.
#' @param level "A", "B", or "C". Used together with keyword.
#' @param children TRUE/FALSE, include children pathways or genes if TRUE.
#' @return data.frame
#' @author ZG Zhao
#' @export
koInfo <- function(kos=NULL, keyword=NULL, level="B", children=FALSE) {
    ss <- TRUE
    if(!is.null(kos)) {
        ss <- ss & keggs$ko %in% kos
    } else if(!is.null(keyword)) {
        s1 <- grepl(keyword, keggs$description, ignore.case = TRUE)
        s2 <- keggs$level %in% level
        ss <- ss & s1 & s2
    }
    if(all(ss)) return(keggs)
    if(children) {
        kos <- keggs$ko[ss]
        ss <- ss | keggs$Parent %in% kos
    }
    dtx <- keggs[ss, ]
    rownames(dtx) <- NULL
    dtx
}

## keggs <- parseKeg("../inst/ko00001.keg")
## str(keggs)
## save(keggs, file="../data/keggs.rda")

## ATHanno <- read.table("~/A.Bioinfo/Data/TAIR10/TAIR10_full_anno.txt", sep="\t",
##                       header=TRUE, row.names=1, stringsAsFactors=FALSE, quote='')
## ATHanno <- ATHanno[, -(c(2, 3, 4, 7))]
## str(ATHanno)
## save(ATHanno, file="../data/ATHanno.rda")
