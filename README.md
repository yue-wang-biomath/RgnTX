# RgnTX

## Instruction
-The Bioconductor package RgnTX allows user-defined restriction of transcriptome areas during the shuffling, and therefore offers high flexibility in the null model to simulate the realistic transcriptome-wide background such as the complex alternative splicing patterns. The setting of null models (randomization) and colocalization measures can be independently chosen from many pre-defined choices when performing statistical colocalization analysis. Importantly, RgnTX also supports the testing of transcriptome elements without clear isoform association, which is often the real scenario due to technical limitations. 

- Function `shiftedZScoreTx` is updated supprting the shifting regions of interest (ROI) over mRNA space (exons). 
- Function `shiftExonTx` is updated for pick regions over mRNA space (exons).
- Files in /inst/extdata/Review are updated for review purpose only.

## 1. Install
To install this package via devtools, please use the following codes.
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yue-wang-biomath/RgnTX", build_vignettes = TRUE)
```

To install this package via BiocManager, please use the following codes.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RgnTX")
```

To view the documentation of RgnTX, please type `browseVignettes("RgnTX")` after installation.

## 2. Basic functions

### - shiftExonTx

```R
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons.tx0 <- exonsBy(txdb)
#Shift regions over the following transcripts
trans.ids <- c("170", "1387", "4113", "10715")
regions <- exons.tx0[trans.ids]
start <- c(3624255, 55158197, 881641, 15694195)
width <- c(200, -200, 200, -200)
shifted_regions <- shiftExonTx(regions, start, width)
```

```
|  FUNCTION shiftExonTx(regions, start, width){
|   FOR EACH region_i in regions
|       Take the start_i in start
|       Take the width_i in width
|       region_target_i <- the region starting from start_i has width_i over region_i
|   ENDFOR
|   RETURN region_target
|  }
```

### - randomizeTx
```R
trans.ids <- c("170", "782", "974", "1364", "1387")
randomResults <- randomizeTx(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             trans_ids = trans.ids, 
                             random_num = 10, 
                             type = "mature", 
                             random_length = 100)
```

## 3. Quick start
### - permTestTx

### - permTestTx_customPick
```R
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
getCDS = function(trans_ids, txdb,...){
  cds.tx0 <- cdsBy(txdb, use.names=FALSE)
  cds.names <- as.character(intersect(names(cds.tx0), trans_ids))
  cds = cds.tx0[cds.names]
  return(cds)
}
RS1 <- randomizeTx(txdb, "all", random_num = 100, random_length = 200, type = "CDS")

permTestTx_results <- permTestTx_customPick(RS1 = RS1,
                                            txdb = txdb,
                                            customPick_function = getCDS,
                                            ntimes = 50,
                                            ev_function_1 = overlapCountsTx,
                                            ev_function_2 = overlapCountsTx)
p1 <- plotPermResults(permTestTx_results, binwidth = 1)
```

### - permTestTxIA_customPick
```R
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
file <- system.file(package="RgnTX", "extdata", "m6A_sites_data.rds")
m6A_sites_data <- readRDS(file)
RS1 <- m6A_sites_data[1:500]

permTestTx_results <- permTestTxIA_customPick(RS1 = RS1,
                                       txdb = txdb,
                                       customPick_function = getStopCodon,
                                       ntimes = 50,
                                       ev_function_1 = overlapCountsTxIA,
                                       ev_function_2 = overlapCountsTx)
p_a <- plotPermResults(permTestTx_results, binwidth = 1)
```

### - permTestTx
```R
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons.tx0 <- exonsBy(txdb)
trans.ids <- sample(names(exons.tx0), 50)

A <- randomizeTx(txdb, trans.ids, random_num = 100, random_length = 100)
B <- c(randomizeTx(txdb, trans.ids, random_num = 50, random_length = 100), 
       A[1:50])

permTestTx_results <- permTestTx(RS1 = A, 
                                 RS2 = B, 
                                 txdb = txdb, 
                                 type = "mature",
                                 ntimes = 50, 
                                 ev_function_1 = overlapCountsTx, 
                                 ev_function_2 = overlapCountsTx)
```

### - plotPermResults

## 4. Shifted z-scores over transcripts

### - shiftedZScoreTx
```R
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
file <- system.file(package="RgnTX", "extdata", "m6A_sites_data.rds")
m6A_sites_data <- readRDS(file)
RS1 <- m6A_sites_data[1:100]
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
permTestTx_results <- permTestTxIA_customPick(RS1 = RS1,
                                              txdb = txdb,
                                              type = "mature",
                                              customPick_function = getStopCodon,
                                              ntimes = 1)
shiftedZScoresTx_results <- shiftedZScoreTx(permTestTx_results,txdb,
                                           type = 'mature',
                                           window = 500,
                                           step = 50,
                                           ev_function_1 = overlapCountsTxIA)
p1 <- plotShiftedZScoreTx(shiftedZScoresTx_results)
```

### - plotShiftedZScoreTx

## 5. Multiple hypothesis tests with Benjamini-Hochberg correction

### - adjustMultipleTesting
```R
file <- system.file(package="RgnTX", "extdata/multi_pvals.rds")
multi_pvals <- readRDS(file)
adjustMultipleTesting(multi_pvals, 0.05)
```
