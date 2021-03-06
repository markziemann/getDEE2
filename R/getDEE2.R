#' getDEE2: Programmatic access to the DEE2 RNA expression dataset
#'
#' Digital Expression Explorer 2 (or DEE2 for short) is a repository of
#' processed RNA-seq data in the form of counts. It was designed so that
#' researchers could undertake re-analysis and meta-analysis of published
#' RNA-seq studies quickly and easily. This package provides an R interface
#' to access these expression data. More information about the DEE2 project
#' can be found at the project homepage (http://dee2.io) and main publication
#' (https://doi.org/10.1093/gigascience/giz022).
#'
#' @docType package
#' @name getDEE2
#' @examples
#' # Example workflow
#' # Fetch metadata
#' mdat <- getDEE2Metadata("celegans")
#' # filter metadata for SRA project SRP009256
#' mdat1 <- mdat[which(mdat$SRP_accession %in% "SRP009256"),]
#' # create a vector of SRA run accessions to fetch
#' SRRvec <- as.vector(mdat1$SRR_accession)
#' # obtain the data as a SummarizedExperiment
#' x <- getDEE2("celegans",SRRvec,metadata=mdat,counts="GeneCounts")
#' # Next, downstream analysis with your favourite Bioconductor tools :)
NULL


getURL <- function(URL, FUN, ..., N.TRIES=1L) {
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
    while (N.TRIES > 0L) {
        result <- tryCatch(FUN(URL, ...), error=identity)
        if (!inherits(result, "error"))
            break
        N.TRIES <- N.TRIES - 1L
    }
    if (N.TRIES == 0L) {
        stop("'getURL()' failed:",
            "\n  URL: ", URL,
            "\n  error: ", conditionMessage(result))
    }
    result
}


#' Get DEE2 Metadata
#'
#' This function fetches the short metadata for the species of interest.
#' @param species A character string matching a species of interest.
#' @param outfile Optional filename.
#' @param ... Additional parameters to be passed to download.file.
#' @return a table of metadata.
#' @keywords metadata
#' @import utils
#' @export
#' @examples
#' ecoli_metadata <- getDEE2Metadata("ecoli")
getDEE2Metadata <- function(species,outfile=NULL, ...){
    orgs = c("athaliana","celegans","dmelanogaster","drerio","ecoli",
    "hsapiens","mmusculus","rnorvegicus","scerevisiae")
    if (species %in% orgs == FALSE ) {
        message(paste("Provided species '",species,"' is not found in the list.
        Check spelling and try again" ,sep=""))
        message(paste("Valid choices are'",paste(orgs,collapse = "', '"),"'."))
    } else {  
        metadataURL = paste("http://dee2.io/metadata/",
        species,"_metadata.tsv.cut",sep="")
        if(is.null(outfile)){
            metadataname = tempfile()
        } else {
            metadataname = outfile
            if (!grepl(".tsv$",metadataname)){
                metadataname = paste0(metadataname,".tsv")
            }
        }
        getURL(URL=metadataURL,FUN=download.file, N.TRIES=1L,
            destfile=metadataname, mode = "wb" , ...)
        mdat <- read.table(metadataname,header=TRUE,quote="",
        sep='\t',fill=TRUE,comment.char = "")
        if (is.null(outfile)){
            unlink(metadataname)
        }
        return(mdat)
    }
}

#' Query Whether a DEE2 Dataset is Available
#'
#' This function sends a query to check whether a dataset is available or not.
#' @param species A character string matching a species of interest.
#' @param SRRvec A character string or vector thereof of SRA run accession
#' numbers.
#' @param metadata optional R object of DEE2 metadata to query.
#' @param ... Additional parameters to be passed to download.file.
#' @return a list of datasets that are present and absent.
#' @keywords query
#' @export
#' @examples
#' x <- queryDEE2("ecoli",c("SRR1067773","SRR5350513"))
queryDEE2 <- function(species, SRRvec,metadata=NULL, ...) {
    if(is.null(metadata)){
        mdat <- getDEE2Metadata(species, ...)
    } else {
        mdat <- metadata
    }
    present <- SRRvec[which(SRRvec %in% mdat$SRR_accession)]
    absent <- SRRvec[-which(SRRvec %in% mdat$SRR_accession)]
    dat <- list("present" = present, "absent" = absent)
    return(dat)
}

#' Load Gene Counts
#'
#' This function loads STAR gene level counts from a downloaded zip file.
#' @param zipname Path to the zipfile.
#' @keywords Load Gene
#' @return a dataframe of gene expression counts. 
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadGeneCounts("mydata.zip")
loadGeneCounts <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("GeneCountMatrix.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,stringsAsFactors=FALSE,
        comment.char = "")
    unlink(TF)
    return(dat)
}

#' Load Transcript Counts 
#'
#' This function loads Kallisto transcript level counts from a downloaded zip
#' file.
#' @param zipname Path to the zipfile.
#' @keywords Load Transcript
#' @return a dataframe of transcript expression counts.
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadTxCounts("mydata.zip")
loadTxCounts <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("TxCountMatrix.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,stringsAsFactors=FALSE,
        comment.char = "")
    unlink(TF)
    return(dat)
}

#' Load Gene Info
#'
#' This function loads gene information. This information includes gene names
#' and lengths which is useful for downstream analysis.
#' @param zipname Path to the zipfile.
#' @keywords Load Gene
#' @return a dataframe of gene information.
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadGeneInfo("mydata.zip")
loadGeneInfo <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("GeneInfo.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,stringsAsFactors=FALSE,
        comment.char = "")
    unlink(TF)
    return(dat)
}

#' Load Transcript Info
#'
#' This function loads transcript information. This information includes
#' transcript lengths, corresponding parent gene accession and gene symbol
#' that might be useful for downstream analysis.
#' @param zipname Path to the zipfile.
#' @keywords Load Transcript
#' @return a dataframe of transcript info
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadTxInfo("mydata.zip")
loadTxInfo <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("TxInfo.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,stringsAsFactors=FALSE,
        comment.char = "")
    unlink(TF)
    return(dat)
}

#' Load Quality Control Info
#'
#' This function loads quality control data. More information about the QC
#' metrics is available from the project github page:
#' https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md
#' @param zipname Path to the zipfile.
#' @keywords Load Qualiy Control QC
#' @return a dataframe of quality control metrics.
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadQcMx("mydata.zip")
loadQcMx <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("QC_Matrix.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,fill=TRUE,
        stringsAsFactors=FALSE,comment.char = "")
    unlink(TF)
    message("For more information about DEE2 QC metrics, visit
    https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md")
    return(dat)
}

#' Load Summary Metadata
#'
#' This function loads the summary metadata, which are the most relevant SRA
#' accession numbers.
#' @param zipname Path to the zipfile.
#' @keywords Load Metadata
#' @return a dataframe of summary metadata.
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadSummaryMeta("mydata.zip")
loadSummaryMeta <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("MetadataSummary.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,quote="",sep='\t',
        fill=FALSE,stringsAsFactors=FALSE,comment.char = "")
    unlink(TF)
    return(dat)
}

#' Load Full Metadata   
#'
#' This function loads the full metadata, which contains many fields.
#' @param zipname Path to the zipfile.
#' @keywords Load Metadata  
#' @return a dataframe of full metadata.      
#' @import utils
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y <- loadFullMeta("mydata.zip")
loadFullMeta <- function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM = FILELIST[grep("MetadataFull.tsv",FILELIST)]
    TF = tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname <- paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,quote="",fill=TRUE,
        sep='\t',stringsAsFactors=FALSE,comment.char = "")
    unlink(TF)
    return(dat)
}

#' Aggregate Transcript Counts to Gene-Level Counts
#'
#' This function converts Kallisto transcript-level expression estimates to
#' gene-level estimates. Counts for each transcript are summed to get an
#' aggregated gene level score.
#' @param x a getDEE2 object.
#' @keywords Aggregate transcript gene
#' @return a dataframe of gene expression counts.
#' @importFrom stats aggregate
#' @export
#' @examples
#' x <- getDEE2("scerevisiae",c("SRR1755149","SRR1755150"),legacy=TRUE)
#' x <- Tx2Gene(x)
Tx2Gene <- function(x){
    y <- merge(x$TxInfo,x$TxCounts,by=0)
    rownames(y) = y$Row.names
    y$Row.names = y$GeneSymbol = y$TxLength = NULL
    yy <- aggregate(. ~ GeneID,y,sum)
    rownames(yy)<-yy$GeneID
    yy$GeneID = NULL
    x <- c(list("Tx2Gene"=yy),x)
}


#' Create summarizedExperiment object
#'
#' This function creates a SummarizedExperiment object from a legacy getDEE2
#' dataset
#' @param x a getDEE2 object.
#' @param counts select "GeneCounts" for STAR based gene counts, "TxCounts" for
#' kallisto transcript level counts or "Tx2Gene" for transcript counts
#' aggregated to gene level. Default is "GeneCounts"
#' @keywords SummarizedExperiment
#' @return a SummarizedExperiment object
#' @import SummarizedExperiment 
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),legacy=TRUE)
#' y <- se(x)
se <- function(x,counts="GeneCounts"){
    if ( (counts == "GeneCounts" ) | is.null(counts) ) {
        myse <- SummarizedExperiment(assays=list(counts=x$GeneCounts),
        colData = data.frame(x$MetadataFull,t(x$QcMx)))
    } else if (counts == "TxCounts") {
        myse <- SummarizedExperiment(assays=list(counts=x$TxCounts),
        colData = data.frame(x$MetadataFull,t(x$QcMx)))
    } else if (counts == "Tx2Gene") {
        if ( length(which(names(x) == "Tx2Gene"))==0 ) {
            x <- Tx2Gene(x)
        }
        myse <- SummarizedExperiment(assays=list(counts=x$Tx2Gene),
        colData = data.frame(x$MetadataFull,t(x$QcMx)))
    } else {
        stop("'Counts' needs to be 'GeneCounts', 'TxCounts' or 'Tx2Gene'")
    }
    myse
}

#' Get DEE2 Gene Expression Data
#'
#' The getDEE2 function fetches gene expression data from the DEE2 database 
#' of RNA sequencing data and returns it as a SummarizedExperiment object.
#' @param species A character string matching the species of interest.   
#' @param SRRvec A character string or vector of SRA run accession numbers.
#' @param counts A string, either 'GeneCounts', 'TxCounts' or 'Tx2Gene'.
#' When 'GeneCounts' is specified, STAR gene level counts are returned.
#' When 'TxCounts' is specified, kallisto transcript counts are returned.
#' When 'Tx2Gene' is specified, kallisto counts aggregated (by sum) on gene
#' are returned. If left blank, "GeneCounts" will be fetched.
#' @param outfile An optional file name for the downloaded dataset.
#' @param metadata (Optional) name of R object for the meta data. Providing
#' the metadata will speed up performance if multiple queries are made in a
#' session. If left blank, the metadata will be fetched once again.
#' @param legacy Whether data should be returned in the legacy (list) format.
#' Default is FALSE. Leave this FALSE if you want to receive data as
#' Summarized experiment.
#' @param baseURL The base URL of the service. Leave this as the default URL
#' unless you want to download from a 3rd party mirror.
#' @param ... Additional parameters to be passed to download.file.
#' @keywords DEE2 RNA-seq database
#' @return a SummarizedExperiment object.
#' @import utils
#' @import SummarizedExperiment 
#' @export
#' @examples
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"))
getDEE2 <- function(species, SRRvec, counts="GeneCounts", metadata=NULL,
    outfile=NULL, legacy=FALSE,
    baseURL="http://dee2.io/cgi-bin/request.sh?", ...){
   SRRvec <- SRRvec[order(SRRvec)]
   if(is.null(metadata)){
        dat1 <- queryDEE2(species, SRRvec)
    } else {
        dat1 <- queryDEE2(species, SRRvec,metadata=metadata)
    }
    absent <- dat1$absent
    present <- dat1$present
    if ( length(present) < 1 ) {
        message("Error. None of the specified SRR accessions are present.")
    } else {
        SRRvec <- gsub(" ","",present)
        llist <- paste0("&x=",paste(SRRvec,collapse = "&x="))
        murl <- paste0(baseURL,"org=",species, llist)
        if(is.null(outfile)){
            zipname = paste(tempfile(),".zip",sep="")
        } else {
            zipname = outfile
            if(!grepl(".zip$",zipname)){
                zipname = paste0(zipname,".zip")
            }
        }
        getURL(URL=murl,FUN=download.file, N.TRIES=1L,
            destfile=zipname, mode = "wb" , ...)
        GeneCounts <- loadGeneCounts(zipname)
        TxCounts <- loadTxCounts(zipname)
        GeneInfo <- loadGeneInfo(zipname)
        TxInfo <- loadTxInfo(zipname)
        QcMx <- loadQcMx(zipname)
        MetadataSummary <- loadSummaryMeta(zipname)
        MetadataFull <- loadFullMeta(zipname)
        dat <- list("GeneCounts" = GeneCounts, "TxCounts" = TxCounts,
        "GeneInfo" = GeneInfo,"TxInfo" = TxInfo , "QcMx" = QcMx,
        "MetadataSummary" = MetadataSummary , "MetadataFull" = MetadataFull ,
        "absent" = absent)
        if(is.null(outfile)){
            unlink(zipname)
        }
        if(length(absent) > 0){
            message(paste0("Warning, datasets not found: '",
            paste(absent,collapse=","),"'"))
        }
        if(legacy == FALSE) {
            if(is.null(counts)) { counts="GeneCounts" }
            dat2 <- se(dat,counts=counts)
            return(dat2)
        }
        if(legacy == TRUE) {
            return(dat)
        }
    }
}


#' Get a table of all completed projects at DEE2
#'
#' This function fetches a table listing all completed projects that
#' are available at DEE2
#' @param species A character string matching a species of interest.
#' @return a table of project bundles available at DEE2.io/huge
#' @keywords metadata
#' @import utils
#' @import htm2txt
#' @export
#' @examples
#' bundles <- list_bundles("celegans")
list_bundles <- function(species){
    orgs = c("athaliana","celegans","dmelanogaster","drerio","ecoli",
    "hsapiens","mmusculus","rnorvegicus","scerevisiae")
    if (species %in% orgs == FALSE ) {
        message(paste("Provided species '",species,"' is not found in the list.
        Check spelling and try again" ,sep=""))
        message(paste("Valid choices are'",paste(orgs,collapse = "', '"),"'."))
    } else {
        URLBASE = paste("http://dee2.io/huge/",species,"/",sep="")
        BUNDLES_FILE = tempfile()
        getURL(URL=URLBASE,FUN=download.file, N.TRIES=1L,
            destfile=BUNDLES_FILE, mode = "wb" )
        bundles <- htm2txt(readLines(BUNDLES_FILE))
        unlink(BUNDLES_FILE)
        bundles <- bundles[grep("RP",bundles)]
        bundles <- t(data.frame(vapply(X=bundles,FUN=function(x) {
            strsplit(x," ")  },list(1))))
        bundles <- as.data.frame(bundles,stringsAsFactors=FALSE)
        rownames(bundles) <- seq(bundles[,1])
        bundles$SRP <- unlist(lapply(strsplit(bundles[,1],"_"),"[[",1))
        bundles$GSE <- unlist(lapply(strsplit(bundles[,1],"_"),"[[",2))
        bundles$GSE <- gsub(".zip","",bundles$GSE)
        colnames(bundles) <- c("file_name","date_added","time_added",
            "file_size","SRP_accession","GSE_accession")
        attributes(bundles)$species <- species
        return(bundles)
    }
}


#' Query whether a project bundle is available from DEE2
#'
#' This function sends a query to check whether a dataset is available or not.
#' @param species A character string matching a species of interest.
#' @param query A character string, such as  the SRA project accession number
#' or the GEO series accession number
#' @param col the column name to be queried, usually "SRP_accession" for SRA
#' project accession or "GSE_accession" for GEO series accession.
#' @param bundles optional table of previously downloaded bundles.
#' @return a list of datasets that are present and absent.
#' @keywords query
#' @export
#' @examples
#' query_bundles("celegans", c("SRP133403","SRP133439"), col = "SRP_accession")
query_bundles <- function(species,query,col,bundles=NULL){
    cols <- c("file_name","date_added","time_added","file_size",
        "SRP_accession","GSE_accession")
    if (col %in% cols == FALSE ) {
        message(paste("Provided column '",col,"' is not found in the list.
        Check spelling and try again" ,sep=""))
        message(paste("Valid choices are'",paste(cols,collapse = "', '"),"'."))
    } else {
        if(is.null(bundles)){
            bundles <- list_bundles(species)
        }
        present <- query[which(query %in% bundles[,col])]
        absent <- query[-which(query %in% bundles[,col])]
        res <- bundles[which(bundles[,col] %in% query),]
        dat <- list("present" = present, "absent" = absent)
        return(dat)
    }
}


#' Get a DEE2 project bundle
#'
#' The getDEE2_bundle function fetches gene expression data from DEE2.
#' This function will only work if all SRA runs have been successfully 
#' processed for an SRA project. This function returns a 
#' SummarizedExperiment object.
#' @param species A character string matching the species of interest.   
#' @param query A character string, such as  the SRA project accession number
#' or the GEO series accession number
#' @param col the column name to be queried, usually "SRP_accession" for SRA
#' project accession or "GSE_accession" for GEO series accession.
#' @param counts A string, either 'GeneCounts', 'TxCounts' or 'Tx2Gene'.
#' When 'GeneCounts' is specified, STAR gene level counts are returned.
#' When 'TxCounts' is specified, kallisto transcript counts are returned.
#' When 'Tx2Gene' is specified, kallisto counts aggregated (by sum) on gene
#' are returned. If left blank, "GeneCounts" will be fetched.
#' @param bundles optional table of previously downloaded bundles.
#' providing this will speed up performance if multiple queries are made in a
#' session. If left blank, the bundle list will be fetched again.
#' @param legacy Whether data should be returned in the legacy (list) format.
#' Default is FALSE. Leave this FALSE if you want to receive data as
#' Summarized experiment.
#' @param baseURL The base URL of the service. Leave this as the default URL
#' unless you want to download from a 3rd party mirror.
#' @param ... Additional parameters to be passed to download.file.
#' @keywords DEE2 RNA-seq database
#' @return a SummarizedExperiment object.
#' @import utils
#' @import SummarizedExperiment 
#' @export
#' @examples
#' x <- getDEE2_bundle("celegans", "SRP133403",col="SRP_accession")
getDEE2_bundle <- function(species, query, col, counts="GeneCounts", 
    bundles=NULL, legacy=FALSE, baseURL="http://dee2.io/huge/", ...){
    if(is.null(bundles)){
        bundles <- list_bundles(species)
    }
    dat1 <- query_bundles(species, query, col, bundles=bundles)
    absent <- dat1$absent
    present <- dat1$present
    if ( length(present) < 1 ) {
        stop("Error. None of the specified accessions are present.")
    }
    if ( length(query) > 1 ) {
        stop("Only one query dataset at a time")
    }
    zipname = paste(tempfile(),".zip",sep="")
    fname = bundles[which(bundles[,col] %in% query),1]
    murl = paste(baseURL,species,"/",fname,sep="")
    getURL(URL=murl,FUN=download.file, N.TRIES=1L,
        destfile=zipname, mode = "wb" , ...)
    GeneCounts <- loadGeneCounts(zipname)
    TxCounts <- loadTxCounts(zipname)
    GeneInfo <- loadGeneInfo(zipname)
    TxInfo <- loadTxInfo(zipname)
    QcMx <- loadQcMx(zipname)
    MetadataSummary <- loadSummaryMeta(zipname)
    MetadataFull <- loadFullMeta(zipname)
    dat <- list("GeneCounts" = GeneCounts, "TxCounts" = TxCounts,
    "GeneInfo" = GeneInfo,"TxInfo" = TxInfo , "QcMx" = QcMx,
    "MetadataSummary" = MetadataSummary , "MetadataFull" = MetadataFull ,
    "absent" = absent)
    unlink(zipname)
    if(length(absent)>0){
        message(paste0("Warning, datasets not found: '",
        paste(absent,collapse=","),"'"))
    }
    if(legacy == FALSE) {
        if(is.null(counts)) { counts="GeneCounts" }
        dat2<-se(dat,counts=counts)
        return(dat2)
    }
    if(legacy == TRUE) {
        return(dat)
    }
}

#' Summarized run data to experiments
#'
#' Sometimes, each SRA experiment data is represented in two or more runs and
#' they need to be aggregated.
#' @param x a getDEE2 object.
#' @param counts select "GeneCounts" for STAR based gene counts, "TxCounts" for
#' kallisto transcript level counts or "Tx2Gene" for transcript counts
#' aggregated to gene level. Default is "GeneCounts"
#' @return a dataframe with gene expression data summarised to SRA experiment
#' accession numbers rather than run accession numbers.
#' @export
#' @examples
#' x <- getDEE2("ecoli",c("SRR1613487","SRR1613488"),legacy=TRUE)
#' y <- srx_agg(x)
srx_agg <- function(x,counts="GeneCounts") {
    if (is(object = x,class2 = "SummarizedExperiment")) {
        stop("This function works only for data obtained with legacy=TRUE")
    }
    IDX=which(names(x) %in% "GeneCounts")
    mds<-x$MetadataSummary
    n=nrow(x[[IDX]])
    SRX_dat <- vapply(X=unique(mds$SRX_accession) ,function(srx) {
        srrs<-rownames(mds)[which(mds$SRX_accession %in% srx)]
        if (length(srrs)>1) {
            rowSums(x[[IDX]][,srrs])
        } else {
            x[[IDX]][,srrs]
        }
    } , numeric(n))
    rownames(SRX_dat) <- rownames(x[[IDX]])
    SRX_dat
}
