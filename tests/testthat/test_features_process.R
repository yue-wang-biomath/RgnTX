library(RgnTX)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Load data the following tests need.
# Load data the following tests need.
file1 <- system.file(package = "RgnTX", "extdata/randomRegionSet1.rds")
randomRegionSet1 <- readRDS(file1)

test_that("Test the GRanges2GRangesList and the GRangesList2GRanges function", {
    randomRegionSet1_gr <- GRangesList2GRanges(randomRegionSet1)
    expect_s4_class(randomRegionSet1_gr, "GRanges")

    randomRegionSet1_list <- GRanges2GRangesList(randomRegionSet1_gr)
    expect_s4_class(randomRegionSet1_list, "GRangesList")
})

test_that("Test the shiftTx function", {
    trans.id.pstv <- c("170", "782", "974", "1364", "1387")
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    cds.tx0 <- cdsBy(txdb, use.names = FALSE)
    cds.p <- cds.tx0[trans.id.pstv]
    width <- 200
    start <- as.numeric(max(end(cds.p)))
    R.cds.last200 <- shiftTx(cds.p, start = start, width = width, direction = "left", strand = "+")
    expect_s4_class(R.cds.last200, "GRanges")

    R.cds.last200.list <- GRanges2GRangesList(R.cds.last200)
    expect_equal(length(R.cds.last200.list), 5)
})
