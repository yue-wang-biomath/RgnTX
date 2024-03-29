# RgnTX

## Introduction
RgnTX is an R/Bioconductor tool for the colocalization analysis of transcriptome elements with permutation tests. Different from existing approaches, RgnTX directly takes advantage of transcriptome annotation, and offers high flexibility in the null model to simulate realistic transcriptome-wide background, such as the complex alternative splicing patterns. The setting of null models (randomization) and colocalization measures can be independently chosen from many pre-defined choices when performing statistical colocalization analysis. Importantly, RgnTX supports the testing of transcriptome elements without clear isoform association, which is often the real scenario due to technical limitations. 

- Function `shiftedZScoreTx` is updated supprting the shifting regions of interest (ROI) over mRNA space (exons). 
- Function `shiftExonTx` is updated for picking regions over mRNA space (exons).
- Files in /inst/extdata/Review are updated for review purpose only (23, Aug 13 updated).
- The ReadMe file is updated providing basic examples and corresponding pseudocodes.

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
The function `shiftExonTx` is able to shift regions over mRNA space, i.e., spanning multiple intervals (exons). It requires three inputs: 'regions', 'start' and 'width'. 

The following example shows `shiftExonTx` pick regions of width 200 over 'regions' starting from each 'start' position. The width value is positive/negative if users would like to move the starting positions to the right (3')/left (5'). 

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

```R
> shifted_regions 
GRangesList object of length 4:
$`170`
GRanges object with 2 ranges and 0 metadata columns:
      seqnames          ranges strand
         <Rle>       <IRanges>  <Rle>
  [1]     chr1 3624255-3624355      +
  [2]     chr1 3638585-3638683      +
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths

$`1387`
GRanges object with 2 ranges and 0 metadata columns:
      seqnames            ranges strand
         <Rle>         <IRanges>  <Rle>
  [1]     chr1 55152023-55152121      +
  [2]     chr1 55158097-55158197      +
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths

$`4113`
GRanges object with 2 ranges and 0 metadata columns:
      seqnames        ranges strand
         <Rle>     <IRanges>  <Rle>
  [1]     chr1 881553-881641      -
  [2]     chr1 880923-881033      -
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths

$`10715`
GRanges object with 4 ranges and 0 metadata columns:
      seqnames            ranges strand
         <Rle>         <IRanges>  <Rle>
  [1]     chr2 15701312-15701351      -
  [2]     chr2 15698704-15698758      -
  [3]     chr2 15696907-15696943      -
  [4]     chr2 15694195-15694262      -
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths
  
 > width(shifted_regions)
IntegerList of length 4
[["170"]] 101 99
[["1387"]] 99 101
[["4113"]] 89 111
[["10715"]] 40 55 37 68
```

```
|  FUNCTION shiftExonTx(regions, start, width){
|   FOR EACH region_i in regions
|       Take the start_i in start
|       Take the width_i in width
|       Extract target_region_i: the region starting from start_i has width_i in region_i
|   ENDFOR
|   Obtain target_region: collection of target_region_i
|
|   RETURN target_region
|  }
```

### - randomizeTx
The function `randomizeTx` randomly picks regions over specified transcripts. 
The following example shows `randomizeTx` picks 10 random regions from a set formed by five mRNAs, with each region being picked having width 100. 

```R
trans.ids <- c("170", "782", "974", "1364", "1387")
randomResults <- randomizeTx(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             trans_ids = trans.ids, 
                             random_num = 10, 
                             type = "mature", 
                             random_length = 100)
```

```
|  FUNCTION randomizeTx(txdb, trans_ids, random_num, random_length, type, N){
|   FOR j in 1 to N
|
|       FOR i in 1 to random_num
|       Take a random element in trans_ids as trans_id_i
|       Take the i-th element in random_length as random_length_i
|
|           IF type is one from "mature"/"full"/"fiveUTR"/"CDS"/"threeUTR"
|               region_i <- Pick corresponding type of region from txdb that has transcript id trans_id_i 
|           ENDIF
|
|       Randomly pick ramdom_region_i: a random region on the transcript region_i having width random_length_i
|       ENDFOR
|       Obtain ramdom_region: collection of ramdom_region_i
|
|   ENDFOR
|   Obtain ramdom_region_list: collection of ramdom_region
|
|   RETURN ramdom_region_list
|  }
```

## 3. Quick start
### - permTestTx
The function `permTestTx` performs permutation test to evaluate association between region set 'RS1' and region set 'RS2'. 

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

```
|  FUNCTION permTestTx(RS1, RS2, txdb, type, ntimes, ev_function_1, ev_function_2){
|
|   observed_overlap: overlapping counts between RS1 and RS2, calculated by ev_function_1
|   
|   FOR i in 1 to LENGTH(RS1)  
|       Take the i-th element in RS1 as RS1_i 
|       random_region_i: a randomized region of RS1_i on the transcript where RS1_i is located. 
|                        The transcript type is specified by the type variable.
|   ENDFOR     
|   random_region: collection of random_region_i
|   random_overlap: overlapping counts between RS and random_region, calculated by ev_function_1
|   
|   Repeat the above procedures ntimes
|   Obtain random_region_list: collection of random_region
|   Obtain random_overlap_list: collection of random_overlap
|
|   # Obtain pvalue and zscore       
|   FOR j in 1 to ntimes
|       Calculate pvalue and zscore based on random_overlap_list and observed_overlap
|   ENDFOR
|
|   RETURN observed_overlap, ramdom_region_list, random_overlap_list, pvalue, zscore
|  }
```

### - permTestTx_customPick
The function `permTestTx_customPick` performs permutation test to evaluate association between features (RS) and transcriptome regions of interest (ROI). The ROI is not directly provided, but generated by the input 'customPick_function'. 

The following example shows `permTestTx_customPick` tests the association between region set 'RS1' with the CDS. 

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
```

```
|  FUNCTION permTestTx_customPick(RS1, txdb, customPick_function, ntimes, ev_function_1, ev_function_2){
|
|   # RS means region set
|   # ROI denotes regions of interest
|   
|   # Obtain ROI 
|   FOR i in 1 to LENGTH(RS1)  
|       Take the i-th element in RS1 as RS1_i 
|       customPick_function takes RS1_i as input and returns ROI_i
|   ENDFOR
|   Obtain ROI: collection of ROI_i
|
|   # Obtain the observed overlap
|   observed_overlap: overlapping counts between RS and ROI, calculated by ev_function_1
|
|   # Obtain multiple random overlaps
|   FOR i in 1 to LENGTH(RS1)  
|       Take the i-th element in RS1 as RS1_i 
|       random_region_i: a randomized region of RS1_i on the transcript where RS1_i is located.
|   ENDFOR
|   random_region: collection of random_region_i
|   random_overlap: overlapping counts between RS and random_region, calculated by ev_function_1
|   
|   Repeat the above procedures ntimes
|   Obtain random_region_list: collection of random_region
|   Obtain random_overlap_list: collection of random_overlap
|
|   # Obtain pvalue and zscore       
|   FOR j in 1 to ntimes
|       Calculate pvalue and zscore based on random_overlap_list and observed_overlap
|   ENDFOR
|
|   RETURN ROI, ramdom_region_list observed_overlap, random_overlap_list, pvalue, zscore
|  }
```

### - permTestTxIA_customPick
The function `permTestTxIA_customPick` performs permutation test to evaluate association between genomic features (RS) that one only know its genomic positions and do not know which specific transcript is should be aligned to, with transcriptome regions of interest (ROI).The ROI is not directly provided, but generated by the input 'customPick_function'. 

The following example shows `permTestTxIA_customPick` tests the association between m6A sites (which specific isoform it is located on is not known) with the stop codons. 

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
```

```
|  FUNCTION permTestTxIA_customPick(RS1, txdb, customPick_function, ntimes, ev_function_1, ev_function_2){
|
|   # RS means region set
|   # ROI denotes regions of interest
|   
|   FOR i in 1 to LENGTH(RS1)  
|       Take the i-th element in RS1 as RS1_i 
|       Find K isoforms that RS1_i can be aligned to
|           FOR k in 1 to K
|               ROI_i_k : customPick_function picks region over the k-th isoform
|               random_region_i_k: pick a random region having the same width with RS1_i over the k-th isoform
|           ENDFOR
|           Obtain ROI_i: collection of ROI_i_k  
|           Obtain random_region_i collection random_region_i_k
|       
|   observed_overlap_i: overlapping counts between RS1_i and ROI_i / K
|   random_overlap_i: overlapping counts between RS1_i and random_region_i_k / K      
|   
|   ENDFOR
|   observed_overlap: collection of observed_overlap_i  
|   random_overlap: collection of random_overlap_i  
|   
|   Repeat the above procedures ntimes
|   Obtain random_region_list: collection of random_region
|   Obtain random_overlap_list: collection of random_overlap
|
|   # Obtain pvalue and zscore       
|   FOR j in 1 to ntimes
|       Calculate pvalue and zscore based on random_overlap_list and observed_overlap
|   ENDFOR
|
|   RETURN ROI, ramdom_region_list observed_overlap, random_overlap_list, pvalue, zscore
|  }
```

### - plotPermResults 
The function `plotPermResults` accepts results from the output of the above permutation test functions and returns a figure visualizing permutation results. 
```R
p_a <- plotPermResults(permTestTx_results, binwidth = 1)
```

<img src = 'https://github.com/yue-wang-biomath/RgnTX/blob/master/vignettes/figures/section5.6.jpg' width = '500px'> 

## 4. Shifted z-scores over transcripts

### - shiftedZScoreTx
Different from the `localZScore` function in regioneR, `shiftZScoreTX` in RgnTX can shift the transcriptome ROIs on mRNA rather than just on the genome coordinates. In the following example, we tested association of m6A and stop codon regions with window 500 and step 50. `shiftZScoreTX` shifts stop codons over the mRNA space (exons). The stop codon region tested here is formed by the last 100 nt CDS and the first 100 nt 3'UTR. As the results show, we can see m6A sites tested here are enriched more on 3'UTR than on CDS, since the peak is deviated to the right (3') direction. 

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
                                              ntimes = 50)
shiftedZScoresTx_results <- shiftedZScoreTx(permTestTx_results,txdb,
                                           type = 'mature',
                                           window = 500,
                                           step = 50,
                                           ev_function_1 = overlapCountsTxIA)
```

### - plotShiftedZScoreTx
The function `plotPermResults` accepts results from the output of the `shiftedZScoreTx` function and returns a figure visualizing shifted z-scores. 

```R
p1 <- plotShiftedZScoreTx(shiftedZScoresTx_results)
```

<img src = 'https://github.com/yue-wang-biomath/RgnTX/blob/master/vignettes/figures/section9.jpg' width = '400px'> 

When more m6A sites and larger window are involved, results are shown as follows.

```R
RS1 <- m6A_sites_data[1:500]
permTestTx_results2 <- permTestTxIA_customPick(RS1 = RS1,
                                              txdb = txdb,
                                              type = "mature",
                                              customPick_function = getStopCodon,
                                              ntimes = 50)
shiftedZScoresTx_results2 <- shiftedZScoreTx(permTestTx_results2,txdb,
                                           type = 'mature',
                                           window = 800,
                                           step = 50,
                                           ev_function_1 = overlapCountsTxIA)
p2 <- plotShiftedZScoreTx(shiftedZScoresTx_results2)
```

<img src = 'https://github.com/yue-wang-biomath/RgnTX/blob/master/vignettes/figures/section9.1.jpg' width = '400px'> 

## 5. Multiple hypothesis tests with Benjamini-Hochberg correction

### - adjustMultipleTesting
This function provides Benjamini-Hochberg correction for adjusting p-values from multiple hypothesis testing. It returns a data frame object with columns of rank, pvals, adjusted_pvals and decision whether to reject null hypothesis. It also reports the proportion of the tests that reject the null hypothesis H0. 
```R
file <- system.file(package="RgnTX", "extdata/multi_pvals.rds")
multi_pvals <- readRDS(file)
adjustMultipleTesting(multi_pvals[, 1], 0.05)

$Table
   rank       pval adjusted_pval reject_H0
1     1 0.00990099    0.01430143       Yes
2     2 0.00990099    0.01430143       Yes
...
$Proportion
[1] 0.8791209
```
