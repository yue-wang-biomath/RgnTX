library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
getreduceSpace = function(regions.trans){
    regions.trans <- unlist(regions.trans)
    regions.trans <-reduce(regions.trans)
    regions.trans.df <- data.frame(regions.trans)
    regions.trans.df.p = regions.trans.df[regions.trans.df[, 'strand'] == '+',]
    regions.trans.df.n = regions.trans.df[regions.trans.df[, 'strand'] == '-',]
    regions.trans.p = regions.trans[regions.trans.df[, 'strand'] == '+']
    regions.trans.n = regions.trans[regions.trans.df[, 'strand'] == '-']
    regions.trans = c(regions.trans.p, regions.trans.n)

    regions.trans.seqnames.p <- unique(regions.trans.df.p[, 'seqnames'])

    regions.group.p <- lapply(regions.trans.df.p[,1], function(x){
        return(which(regions.trans.seqnames.p == x))
    })

    group.p.total <- length(regions.trans.seqnames.p)

    regions.trans.seqnames.n <- unique(regions.trans.df.n[, 'seqnames'])
    regions.group.n <- lapply(regions.trans.df.n[,1], function(x){
        return(which(regions.trans.seqnames.n == x) + group.p.total)
    })

    regions.group.p <- unlist(regions.group.p)
    regions.group.n <- unlist(regions.group.n)

    regions.group <- c(regions.group.p, regions.group.n)
    regions.trans$group <- regions.group
    regions.trans <- GRanges2GRangesList(regions.trans)

    names(regions.trans) <- as.character(unique(regions.group))
    return(regions.trans)}
    regions.exons <- exonsBy(txdb)
    n_perm <- 100
    N <- 100
    # case 1
    m6A_sites_data <- readRDS("~/RgnTX/inst/extdata/m6A_sites_data.rds")
    regions1.A = m6A_sites_data[1:100]
    regions1.D <- randomizeFeaturesTxIA(regions1.A, txdb, N = n_perm)
    trans.info <- getTransInfo(regions1.A, txdb)
    trans.id <- trans.info[, 'trans_ID']
    regions1.C <- getStopCodon(trans.id, txdb)
    orig.ev1 <- overlapCountsTxIA(regions1.A, regions1.C, txdb)
    random.ev1 <- lapply(regions1.D, function(x){
        return(overlapCountsTx(x, regions1.C))})
    random.ev1 <- unlist(random.ev1)

    # case 2
    #regions.exons <- exonsBy(txdb)
    #trans.id <- names(regions.exons)
    #regions2.C <- getStopCodon(trans.id, txdb)
    orig.ev2 <- overlapCountsTxIA(regions1.A, regions2.C, txdb)
    trans.id <- names(regions.exons)
    regions2.D = list()

    for (i in 1:N) {
        trans.id.sample <- sample(trans.id, length(regions1.A), replace = TRUE)
        regions2.D.sample <- regions.exons[(trans.id.sample)]
        regions2.D[i] <- randomizeTransByOrder(regions2.D.sample, 2)
    }
    random.ev2 <- lapply(regions2.D, function(x){
        return(overlapCountsTx(x, regions2.C))})
    random.ev2 <- unlist(random.ev2)

    # case 3
    trans.info <- getTransInfo(regions1.A, txdb)
    trans.id <- trans.info[, 'trans_ID']
    regions3.D <- regions.exons[unique(trans.id)]
    regions3.D.reduced  <- getreduceSpace(regions3.D)
    regions3.F = list()

    width.regions3.D <- sum(width(regions3.D.reduced))
    board <- cumsum(width.regions3.D/sum(width.regions3.D))

    for (i in 1:N) {
        index.runif <- runif(length(regions1.A))
        index.sample <- lapply(index.runif, function(x){return(min(which(x<board)))})
        index.sample <-  unlist(index.sample)
        regions3.D.reduced.sample <- regions3.D.reduced[index.sample]
        regions3.F[i] <- randomizeTransByOrder(regions3.D.reduced.sample, 2)
    }

    regions3.C <- getreduceSpace(regions1.C)
    regions3.C <- unlist(regions3.C)
    orig.ev3 <- numOverlaps(regions1.A, regions3.C)

    random.ev3 <- lapply(regions3.F, function(x){
        return(numOverlaps(x, regions3.C))})
    random.ev3 <- unlist(random.ev3)

    # case 4
    regions4.C <- getreduceSpace(regions2.C)
    regions4.C <- unlist(regions4.C)
    orig.ev4 <- numOverlaps(regions1.A, regions4.C)
    trans.id <- names(regions.exons)
    regions4.F = list()


    for (i in 1:N) {
        trans.id.sample <- sample(trans.id, length(regions1.A), replace = TRUE)
        regions4.D <- regions.exons[unique(trans.id.sample)]
        regions4.D.reduced  <- getreduceSpace(regions4.D)

        width.regions4.D <- sum(width(regions4.D.reduced))
        board <- cumsum(width.regions4.D/sum(width.regions4.D))
        index.runif <- runif(length(regions1.A))
        index.sample <- lapply(index.runif, function(x){return(min(which(x<board)))})
        index.sample <-  unlist(index.sample)

        regions4.D.reduced.sample <-  regions4.D.reduced[index.sample]
        regions4.F[i] <- randomizeTransByOrder(regions4.D.reduced.sample, 2)
    }

    random.ev4 <- lapply(regions4.F, function(x){
        return(numOverlaps(x, (regions4.C)))})
    random.ev4 <- unlist(random.ev4)

    pvals <- c(
        (sum(orig.ev1 <= random.ev1, na.rm = TRUE) + 1) / (N + 1),
        (sum(orig.ev2 <= random.ev2, na.rm = TRUE) + 1) / (N + 1),
        (sum(orig.ev3 <= random.ev3, na.rm = TRUE) + 1) / (N + 1),
        (sum(orig.ev4 <= random.ev4, na.rm = TRUE) + 1) / (N + 1))
