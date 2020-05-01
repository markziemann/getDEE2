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
#' mdat<-getDee2Metadata("celegans")
#' # filter metadata for SRA project SRP009256
#' mdat1<-mdat[which(mdat$SRP_accession %in% "SRP009256"),]
#' # create a vector of SRA run accessions to fetch
#' SRRvec<-as.vector(mdat1$SRR_accession)
#' # obtain the data as a SummarizedExperiment
#' x<-getDEE2("celegans",SRRvec,metadata=mdat,counts="GeneCounts")
#' # Next, downstream analysis with your favourite Bioconductor tools :)
NULL

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
#' ecoli_metadata<-getDee2Metadata("ecoli")
getDee2Metadata<-function(species,outfile=NULL, ...){
    orgs=c("athaliana","celegans","dmelanogaster","drerio","ecoli","hsapiens",
    "mmusculus","rnorvegicus","scerevisiae")
    if (species %in% orgs == FALSE ) {
        message(paste("Provided species '",species,"' is not found in the list. 
        Check spelling and try again" ,sep=""))
        message(paste("Valid choices are'",paste(orgs,collapse = "', '"),"'."))
    } else {  
        metadataURL=paste("http://dee2.io/metadata/",
        species,"_metadata.tsv.cut",sep="")
        if(is.null(outfile)){
            metadataname=tempfile()
        } else {
            metadataname=outfile
            if (!grepl(".tsv$",metadataname)){
                metadataname=paste0(metadataname,".tsv")
            }
        }
        download.file(metadataURL, destfile=metadataname, mode = "wb", ...)
        mdat<-read.table(metadataname,header=TRUE,quote="",
        sep='\t',fill=TRUE)
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
#' x<-queryDee2("ecoli",c("SRR1067773","SRR5350513"))
queryDee2<-function(species, SRRvec,metadata=NULL, ...) {
    if(is.null(metadata)){
        mdat<-getDee2Metadata(species, ...)
    } else {
        mdat<-metadata
    }
    present<-SRRvec[which(SRRvec %in% mdat$SRR_accession)]
    absent<-SRRvec[-which(SRRvec %in% mdat$SRR_accession)]
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadGeneCounts("mydata.zip")
loadGeneCounts<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("GeneCountMatrix.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE)
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadTxCounts("mydata.zip")
loadTxCounts<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("TxCountMatrix.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE)
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadGeneInfo("mydata.zip")
loadGeneInfo<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("GeneInfo.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE)
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadTxInfo("mydata.zip")
loadTxInfo<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("TxInfo.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE)
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadQcMx("mydata.zip")
loadQcMx<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("QC_Matrix.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,fill=TRUE)
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadQcMx("mydata.zip")
loadSummaryMeta<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("MetadataSummary.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,quote="",sep='\t',fill=FALSE)
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"),outfile="mydata.zip")
#' y<-loadQcMx("mydata.zip")
loadFullMeta<-function(zipname){
    FILELIST <- unzip(zipname,list=TRUE)$Name
    CM=FILELIST[grep("MetadataFull.tsv",FILELIST)]
    TF=tempfile()
    unzip(zipname, files = CM, exdir = tempdir() )
    mxname<-paste0(tempdir(),"/",CM)
    file.rename(mxname,TF)
    dat <- read.table(TF,row.names=1,header=TRUE,quote="",fill=TRUE,sep='\t')
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
#' x<-getDEE2_legacy("scerevisiae",c("SRR1755149","SRR1755150"))
#' x<-Tx2Gene(x)
Tx2Gene<-function(x){
    y<-merge(x$TxInfo,x$TxCounts,by=0)
    rownames(y)=y$Row.names
    y$Row.names=y$GeneSymbol=y$TxLength=NULL
    yy<-aggregate(. ~ GeneID,y,sum)
    rownames(yy)<-yy$GeneID
    yy$GeneID=NULL
    x<-c(list("Tx2Gene"=yy),x)
}

#' Get DEE2 Gene Expression Data (Legacy)
#'
#' This function fetches gene expression data from the DEE2 database of RNA
#' sequencing data and returns it as a getDEE2 object.
#' @param species A character string matching the species of interest.   
#' @param SRRvec A character string or vector of SRA run accession numbers.
#' @param outfile An optional file name for the downloaded dataset.
#' @param metadata (Optional) name of R object for the meta data. Providing
#' the metadata will speed up performance if multiple queries are made in a
#' session.
#' @param baseURL The base URL of the service. Leave this as the default URL
#' unless you want to download from a 3rd party mirror.
#' @param ... Additional parameters to be passed to download.file.
#' @keywords DEE2 RNA-seq database
#' @return a getDEE2 object.
#' @import utils
#' @export
#' @examples
#' x<-getDEE2_legacy("ecoli",c("SRR1613487","SRR1613488"))
getDEE2_legacy<-function(species, SRRvec, outfile=NULL, metadata=NULL,
baseURL="http://dee2.io/cgi-bin/request.sh?", ...){
    if(is.null(metadata)){
        dat1<-queryDee2(species, SRRvec)
    } else {
        dat1<-queryDee2(species, SRRvec,metadata=metadata)
    }
    absent<-dat1$absent
    present<-dat1$present
    if ( length(present) < 1 ) {
        stop("Error. None of the specified SRR accessions are present.")
    } 
    SRRvec<-gsub(" ","",present)
    llist<-paste0("&x=",paste(SRRvec,collapse = "&x="))
    murl <- paste0(baseURL,"org=",species, llist)
    if(is.null(outfile)){
        zipname=tempfile()
    } else {
        zipname=outfile
        if(!grepl(".zip$",zipname)){
            zipname=paste0(zipname,".zip")
        }
    }
    download.file(murl, destfile=zipname, mode = "wb", ...)
    GeneCounts<-loadGeneCounts(zipname)
    TxCounts<-loadTxCounts(zipname)
    GeneInfo<-loadGeneInfo(zipname)
    TxInfo<-loadTxInfo(zipname)
    QcMx<-loadQcMx(zipname)
    MetadataSummary<-loadSummaryMeta(zipname)
    MetadataFull<-loadFullMeta(zipname)
    dat <- list("GeneCounts" = GeneCounts, "TxCounts" = TxCounts,
    "GeneInfo" = GeneInfo,"TxInfo" = TxInfo , "QcMx" = QcMx, 
    "MetadataSummary" = MetadataSummary , "MetadataFull" = MetadataFull ,
    "absent" = absent)

    if(is.null(outfile)){
        unlink(zipname)
    }
    if(length(absent)>0){
        message(paste0("Warning, datasets not found: '",
        paste(absent,collapse=","),"'"))
    }
    return(dat)
}

#' Create summarizedExperiment object
#'
#' This function creates a SummarizedExperiment object from a getDEE2 dataset
#' @param x a getDEE2 object.
#' @param counts select "GeneCounts" for STAR based gene counts, "TxCounts" for
#' kallisto transcript level counts or "Tx2Gene" for transcript counts
#' aggregated to gene level. Default is "GeneCounts"
#' @keywords SummarizedExperiment
#' @return a SummarizedExperiment object
#' @import SummarizedExperiment 
#' @export
#' @examples
#' x<-getDEE2_legacy("ecoli",c("SRR1613487","SRR1613488"))
#' y<-se(x)
se<-function(x,counts="GeneCounts"){
    if ( (counts == "GeneCounts" ) | is.null(counts) ) {
        myse<-SummarizedExperiment(assays=list(counts=x$GeneCounts),
        colData=data.frame(x$MetadataFull,t(x$QcMx)))
    } else if (counts == "TxCounts") {
        myse<-SummarizedExperiment(assays=list(counts=x$TxCounts),
        colData=data.frame(x$MetadataFull,t(x$QcMx)))
    } else if (counts == "Tx2Gene") {
        if ( length(which(names(x) == "Tx2Gene"))==0 ) {
            x<-Tx2Gene(x)
        }
        myse<-SummarizedExperiment(assays=list(counts=x$Tx2Gene),
        colData=data.frame(x$MetadataFull,t(x$QcMx)))
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
getDEE2<-function(species, SRRvec, counts="GeneCounts", metadata=NULL,
    outfile=NULL, baseURL="http://dee2.io/cgi-bin/request.sh?", ...){
    if(is.null(metadata)){
        dat1<-queryDee2(species, SRRvec)
    } else {
        dat1<-queryDee2(species, SRRvec,metadata=metadata)
    }
    absent<-dat1$absent
    present<-dat1$present
    if ( length(present) < 1 ) {
        message("Error. None of the specified SRR accessions are present.")
    } else {
        SRRvec<-gsub(" ","",present)
        llist<-paste0("&x=",paste(SRRvec,collapse = "&x="))
        murl <- paste0(baseURL,"org=",species, llist)
        if(is.null(outfile)){
            zipname = paste(tempfile(),".zip",sep="")
        } else {
            zipname=outfile
            if(!grepl(".zip$",zipname)){
                zipname=paste0(zipname,".zip")
            }
        }
        download.file(murl, destfile=zipname, mode = "wb", ...)
        GeneCounts<-loadGeneCounts(zipname)
        TxCounts<-loadTxCounts(zipname)
        GeneInfo<-loadGeneInfo(zipname)
        TxInfo<-loadTxInfo(zipname)
        QcMx<-loadQcMx(zipname)
        MetadataSummary<-loadSummaryMeta(zipname)
        MetadataFull<-loadFullMeta(zipname)
        dat <- list("GeneCounts" = GeneCounts, "TxCounts" = TxCounts,
        "GeneInfo" = GeneInfo,"TxInfo" = TxInfo , "QcMx" = QcMx,
        "MetadataSummary" = MetadataSummary , "MetadataFull" = MetadataFull ,
        "absent" = absent)
        if(is.null(outfile)){
            unlink(zipname)
        }
        if(length(absent)>0){
            message(paste0("Warning, datasets not found: '",
            paste(absent,collapse=","),"'"))
        }
        if(is.null(counts)) { counts="GeneCounts" }
        dat2<-se(dat,counts=counts)
        return(dat2)
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
    orgs=c("athaliana","celegans","dmelanogaster","drerio","ecoli","hsapiens",
    "mmusculus","rnorvegicus","scerevisiae")
    if (species %in% orgs == FALSE ) {
        message(paste("Provided species '",species,"' is not found in the list.
        Check spelling and try again" ,sep=""))
        message(paste("Valid choices are'",paste(orgs,collapse = "', '"),"'."))
    } else {
        URLBASE=paste("http://dee2.io/huge/",species,"/",sep="")
        BUNDLES_FILE=tempfile()
        download.file(URLBASE,destfile=BUNDLES_FILE)
        bundles <- htm2txt(readLines(BUNDLES_FILE))
        unlink(BUNDLES_FILE)
        bundles <- bundles[grep("RP",bundles)]
        bundles <- t(data.frame(sapply(bundles,function(x) { strsplit(x," ")  })))
        bundles <- as.data.frame(bundles,stringsAsFactors=FALSE)
        rownames(bundles) <- seq(bundles[,1])
        bundles$SRP <- sapply(strsplit(bundles[,1],"_"),"[[",1)
        bundles$GSE <- sapply(strsplit(bundles[,1],"_"),"[[",2)
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
            bundles<- list_bundles(species)
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
bundles=NULL, baseURL="http://dee2.io/huge/", ...){
    if(is.null(bundles)){
        bundles <- list_bundles(species)
    }
    dat1<-query_bundles(species, query, col, bundles=bundles)
    absent<-dat1$absent
    present<-dat1$present
    if ( length(present) < 1 ) {
        stop("Error. None of the specified accessions are present.")
    }
    if ( length(query) > 1 ) {
        stop("Only one query dataset at a time")
    }
    zipname = paste(tempfile(),".zip",sep="")
    fname = bundles[which(bundles[,col] %in% query),1]
    murl = paste(baseURL,species,"/",fname,sep="")
    download.file(murl, destfile = zipname, mode = "wb", ...)
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
    if(is.null(counts)) { counts="GeneCounts" }
    dat2<-se(dat,counts=counts)
    return(dat2)
}

