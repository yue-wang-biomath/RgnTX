library(RgnTX)

# Load data the following tests need.
file1 <- system.file(package = "RgnTX", "extdata/randomRegionSet1.rds")
randomRegionSet1 <- readRDS(file1)
file2 <- system.file(package = "RgnTX", "extdata/randomRegionSet2.rds")
randomRegionSet2 <- readRDS(file2)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_that("Test if overlapCountsTx returns a correct value", {
    expect_equal(overlapCountsTx(randomRegionSet1, randomRegionSet2, count_once = TRUE, over_trans = TRUE), 86)
    expect_equal(overlapCountsTx(randomRegionSet1, randomRegionSet2, count_once = FALSE, over_trans = FALSE), 206)
})

test_that("Test if distanceTx returns a correct value", {
    expect_equal(distanceTx(randomRegionSet1, randomRegionSet2), 0)}
    )

test_that("Test if overlapWidthTx returns a correct value", {
    expect_equal(overlapWidthTx(randomRegionSet1, randomRegionSet2), 11183)
})
