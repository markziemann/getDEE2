library("getDEE2")
library("SummarizedExperiment")
library("testthat")

# E. coli 
x<-getDEE2_legacy("ecoli",c("SRR1613487","SRR1613488"))

test_that("eco works", {
    expect_equal( sum(x$GeneCounts) , 20624168 )
})


# A. thaliana bundle
x <- getDEE2_bundle("athaliana", "SRP058781",col="SRP_accession")

test_that("ath bundles work", {
    expect_equal( nrow(assays(x)[[1]]) , 32833 )
    expect_equal( ncol(assays(x)[[1]]) , 32 )
})


# check absent present
dat <- query_bundles("drerio",c("SRP131781","SRP055996","SRXXX"),col="SRP_accession")

test_that("dre bundle query", {
    expect_equal( length(dat$absent) , 1 )
    expect_equal( length(dat$present) , 2 )
})
