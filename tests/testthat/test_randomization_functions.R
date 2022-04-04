library(RgnTX)

# Load data the following tests need.
file1 <- system.file(package = "RgnTX", "extdata/randomRegionSet1.rds")
randomRegionSet1 <- readRDS(file1)
file3 <- system.file(package = "RgnTX", "extdata/m6A_sites_data.rds")
m6A_sites_data <- readRDS(file3)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_that("Test the randomizeTx function", {
    trans.ids <- c("170", "782", "974", "1364", "1387")
    RS1 <- randomizeTx(txdb, trans.ids, random_num = 10, random_length = 100)
    expect_s4_class(RS1, "GRangesList")
    expect_equal(length(RS1), 10)
})

test_that("Test the randomizeFeaturesTx function", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    RS2 <- randomizeFeaturesTx(randomRegionSet1[1:10], txdb, N = 1)
    expect_s4_class(RS2, "GRangesList")
    expect_equal(length(RS2), 10)
})

test_that("Test the randomizeFeaturesTxIA function", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    RS3 <- randomizeFeaturesTxIA(m6A_sites_data[1:10], txdb, N = 1)
    expect_s4_class(RS3, "GRangesList")
    expect_equal(length(RS3), 10)
    RS3.list <- randomizeFeaturesTxIA(m6A_sites_data[1:10], txdb, N = 2)
    expect_equal(length(RS3.list), 2)
    expect_type(RS3.list, "list")
    expect_s4_class(RS3.list[[1]], "GRangesList")
})

test_that("Test the randomizeTransByOrder function", {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    RS4 <- randomizeTransByOrder(randomRegionSet1[1:10], 20)
    expect_s4_class(RS4, "GRangesList")
    expect_equal(length(RS4), 10)
})
