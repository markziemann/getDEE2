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
    CM="GeneCountMatrix.tsv"
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
    CM="TxCountMatrix.tsv"
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
    CM="GeneInfo.tsv"
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
    CM="TxInfo.tsv"
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
    CM="QC_Matrix.tsv"
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
    CM="MetadataSummary.tsv"
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
    CM="MetadataFull.tsv"
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
#' x<-getDEE2("celegans",c("SRR3657716","SRR3657717"))
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
#' @param SRRvec A character string or vector of SRA run accession numbers
#' @param outfile An optional file name for the downloaded dataset.
#' @param metadata An optional file name for the meta data
#' @param baseURL The base URL of the service. Leave this as the default URL
#' unless you want to download from a 3rd party mirror.
#' @param ... Additional parameters to be passed to download.file.
#' @keywords DEE2 RNA-seq database
#' @return a getDEE2 object.
#' @import utils
#' @export
#' @examples
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"))
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
        message("Error. None of the specified SRR accessions are present.")
    } else {
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
#' x<-getDEE2("ecoli",c("SRR1613487","SRR1613488"))
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
#' This function fetches gene expression data from the DEE2 database of RNA
#' sequencing data and returns it as a SummarizedExperiment object.
#' @param species A character string matching the species of interest.   
#' @param SRRvec A character string or vector of SRA run accession numbers
#' @param counts A string, either 'GeneCounts', 'TxCounts' or 'Tx2Gene'.
#' When 'GeneCounts' is specified, STAR gene level counts are returned.
#' When 'TxCounts' is specified, kallisto transcript counts are returned.
#' When 'Tx2Gene' is specified, kallisto counts aggregated (sum)on gene are
#' returned.
#' @param outfile An optional file name for the downloaded dataset.
#' @param metadata An optional file name for the meta data
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
getDEE2<-function(species, SRRvec, outfile=NULL, metadata=NULL, counts=NULL,
baseURL="http://dee2.io/cgi-bin/request.sh?", ...){
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
        if(is.null(counts)) { counts="GeneCounts" }
        dat2<-se(dat,counts=counts)
        return(dat2)
    }
}

