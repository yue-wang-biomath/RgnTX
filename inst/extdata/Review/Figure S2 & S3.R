# Figure S2
library(RgnTX)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
load("~/RgnTX/m7G_methyl.rda")

get5UTR = function(trans_ids, txdb,...){
    tx0 <- fiveUTRsByTranscript(txdb, use.names=FALSE)
    cds.names <- as.character(intersect(names(tx0), trans_ids))
    fiveUTR <- tx0[names]
    return(fiveUTR)
}

permTestTx_results <- permTestTxIA_customPick(RS1 = sample(m7G_methyl, 100),
                                              txdb = txdb,
                                              customPick_function = get5UTR,
                                              ntimes = 100)
# Figure S3
library(RgnTX)
file <- system.file(package="RgnTX", "extdata", "m6A_sites_data.rds")
m6A_sites_data <- readRDS(file)
library(GenomicFeatures)
m6A_chr4 = m6A_sites_data[(seqnames(m6A_sites_data)) == 'chr4']
for(i in 1:50){
    permTestTx_results <- permTestTxIA_customPick(RS1 = sample(m6A_chr4, 100),
                                                  txdb = txdb,
                                                  customPick_function = getStopCodon,
                                                  ntimes = 100)
    saveRDS(permTestTx_results,paste0('~/RgnTX/chr4/permTestTx_results_',i,'.rds'))
}
pval_list <- list()
for(i in 1:50){
    permTestTx_results <- readRDS(paste0('~/RgnTX/chr4/permTestTx_results_',i,'.rds'))
    p_i <- permTestTx_results[[6]]
    pval_list[[i]] <- p_i
}
pval_list <- unlist(pval_list)

if(interactive()){
    hist(pval_list,ylim = range(0, 50))
    topptx(filename = '~/RgnTX/chr4/chr19-pval.pptx')}

if(interactive()){
    hist(pval_adjust_list,ylim = range(0, 50))
    topptx(filename = '~/RgnTX/chr4/chr19-pval-adjust.pptx')}

