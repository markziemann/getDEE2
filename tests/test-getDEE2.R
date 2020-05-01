library("getDEE2")
library("SummarizedExperiment")
library("testthat")

test_that("multiplication works", {
    expect_equal(2 * 2, 4)
})


# E. coli 
x<-getDEE2_legacy("ecoli",c("SRR1613487","SRR1613488"))

test_that("eco works", {
    expect_equal( sum(x$GeneCounts) , 20624168 )
})


# C. elegans
x<-getDEE2_legacy("celegans",c("SRR051935","SRR051934"))

test_that("cel works", {
    expect_equal( sum(x$TxCounts) , 29067 )
})


# M. musculus
x<-getDEE2_legacy("mmusculus",c("SRR1022283","SRR1022284"))

test_that("mmu works", {
    expect_equal( sum(x$GeneCounts) , 14378076 )
    expect_equal( sum(x$TxCounts) , 2872578 )
})


# A. thaliana bundle
x <- getDEE2_bundle("athaliana", "SRP058781",col="SRP_accession")

test_that("ath bundles work", {
    expect_equal( nrow(assays(x)[[1]]) , 32833 )
    expect_equal( ncol(assays(x)[[1]]) , 32 )
})


# D. rerio bundle
# SRP124609_GSE106677.zip
x <- getDEE2_bundle("drerio", "GSE106677",col="GSE_accession")

test_that("dre bundles work", {
    expect_equal( nrow(assays(x)[[1]]) , 35117 )
    expect_equal( ncol(assays(x)[[1]]) , 67 )
})


# check absent present
dat <- query_bundles("drerio",c("SRP131781","SRP055996","SRXXX"),col="SRP_accession")

test_that("dre bundle query", {
    expect_equal( length(dat$absent) , 1 )
    expect_equal( length(dat$present) , 2 )
})
