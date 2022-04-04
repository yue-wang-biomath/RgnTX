library(RgnTX)

# Load data the following tests need.
file1 <- system.file(package = "RgnTX", "extdata/randomRegionSet1.rds")
randomRegionSet1 <- readRDS(file1)
file2 <- system.file(package = "RgnTX", "extdata/randomRegionSet2.rds")
randomRegionSet2 <- readRDS(file2)
file3 <- system.file(package = "RgnTX", "extdata/m6A_sites_data.rds")
m6A_sites_data <- readRDS(file3)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_that("Test if overlapCountsTx returns a correct value", {
    expect_equal(overlapCountsTx(randomRegionSet1, randomRegionSet2, count_once = TRUE, over_trans = TRUE), 65)
    expect_equal(overlapCountsTx(randomRegionSet1, randomRegionSet2, count_once = TRUE, over_trans = FALSE), 166)
    expect_equal(overlapCountsTx(randomRegionSet1, randomRegionSet2, count_once = FALSE, over_trans = TRUE), 72)
    expect_equal(overlapCountsTx(randomRegionSet1, randomRegionSet2, count_once = FALSE, over_trans = FALSE), 256)
})

test_that("Test if distanceTx returns a correct value", {
    expect_equal(distanceTx(randomRegionSet1, randomRegionSet2), 1993.42197, tolerance = 10e-6)
    expect_equal(distanceTx(randomRegionSet1, randomRegionSet2, beta = 0), 1201206.4, tolerance = 10e-6)
})

test_that("Test if overlapWidthTx returns a correct value", {
    expect_equal(overlapWidthTx(randomRegionSet1, randomRegionSet2), 7161)
})
