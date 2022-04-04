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

permTestTx_results <- permTestTx(RS1 = randomRegionSet1, RS2 = randomRegionSet2, txdb = txdb, ntimes = 2)
test_that("Test the permTestTx function", {
    expect_type(permTestTx_results, "list")
    expect_type(permTestTx_results$RSL, "list")
    expect_s4_class(permTestTx_results$RSL[[1]], "GRangesList")
    expect_equal(length(permTestTx_results), 8)
    expect_true(is(permTestTx_results, "permTestTx.results"))
    expect_error(plotPermResults(""), "Argument permTestTx_results must be a permTestTx.results object.")
})

shiftedZScoreTx_results <- shiftedZScoreTx(permTestTx_results, txdb)
test_that("Test the shiftedZScoreTx function", {
    expect_error(shiftedZScoreTx("", txdb), "Argument permTestTx_results must be a permTestTx.results object.")
    expect_type(shiftedZScoreTx_results, "list")
    expect_type(shiftedZScoreTx_results$shifted.z.scores, "double")
    expect_true(is(shiftedZScoreTx_results, "shitedZScoreTx.results"))
    expect_equal(length(shiftedZScoreTx_results), 4)
    expect_error(plotShiftedZScoreTx(""), "Argument shitedZScoresTx_results must be a shitedZScoreTx.results object.")
})

permTestTx_results <- permTestTxIA(RS1 = m6A_sites_data[1:10], RS2 = randomRegionSet2, txdb = txdb, ntimes = 2)
test_that("Test the permTestTxIA_customPick function", {
    expect_type(permTestTx_results, "list")
    expect_type(permTestTx_results$RSL, "list")
    expect_s4_class(permTestTx_results$RSL[[1]], "GRangesList")
    expect_equal(length(permTestTx_results), 8)
    expect_true(is(permTestTx_results, "permTestTx.results"))
    expect_error(plotPermResults(""), "Argument permTestTx_results must be a permTestTx.results object.")
})

getCDS <- function(trans_ids, ...) {
    cds.tx0 <- cdsBy(txdb, use.names = FALSE)
    cds.names <- as.character(intersect(names(cds.tx0), trans_ids))
    cds <- cds.tx0[cds.names]
    return(cds)
}
permTestTx_results <- permTestTx_customPick(randomRegionSet1, txdb = txdb, customPick_function = getCDS, ntimes = 2)
test_that("Test the permTestTx_customPick function", {
    expect_type(permTestTx_results, "list")
    expect_type(permTestTx_results$RSL, "list")
    expect_s4_class(permTestTx_results$RSL[[1]], "GRangesList")
    expect_equal(length(permTestTx_results), 8)
    expect_true(is(permTestTx_results, "permTestTx.results"))
    expect_error(plotPermResults(""), "Argument permTestTx_results must be a permTestTx.results object.")
})

trans.ids <- c("170", "782", "974", "1364", "1387")
RSL <- randomizeTx(txdb, trans.ids, random_num = 20, random_length = 100, N = 2)
permTestTx_results <- permTestTx_customAll(RSL = RSL, RS1 = randomRegionSet1, RS2 = randomRegionSet2)
test_that("Test the permTestTx_customAll function", {
    expect_type(permTestTx_results, "list")
    expect_s4_class(permTestTx_results$RSL[[1]], "GRangesList")
    expect_equal(length(permTestTx_results), 8)
    expect_true(is(permTestTx_results, "permTestTx.results"))
    expect_error(plotPermResults(""), "Argument permTestTx_results must be a permTestTx.results object.")
})
