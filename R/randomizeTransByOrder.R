#' Randomize features into transcriptome
#'
#' @export randomizeTransByOrder
#' @importFrom stats runif
#'
#' @description This function receives a \code{GRangesList} object and picks a random region within each list element of this object. The length of the region to be picked is decided by the input \code{random_length} argument.
#' @usage randomizeTransByOrder(regions_A, random_length = 20)
#'
#' @param regions_A A \code{GRangesList} object. The name of each list element should be the corresponding transcript id.
#' @param random_length A \code{numeric} object.
#'
#' @return A \code{GRangesList} object. The name of each list element should be the corresponding transcript id.
#'
#' @seealso \code{\link{randomizeTx}}, \code{\link{randomizeFeaturesTx}}, \code{\link{randomizeFeaturesTxIA}}
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' exons.tx0 <- exonsBy(txdb)
#' trans.ids <- sample(names(exons.tx0), 500)
#' regions.A <- exons.tx0[trans.ids]
#' RS <- randomizeTransByOrder(regions.A, random_length = 20)
randomizeTransByOrder <-
    function(regions_A,
             random_length = 20){
        random.length <- random_length
        regions.A <- regions_A
        A.widths <- width(regions.A)
        A.cumsum <- cumsum(A.widths)
        A.length <- max(A.cumsum)
        random.num <- length(regions.A)

        p.index <- which(data.frame(unique(strand(regions.A)))[, 'value'] == '+')
        n.index <- which(data.frame(unique(strand(regions.A)))[, 'value'] == '-')
        suppressWarnings(
            dist <- round(runif(random.num)*(A.length - random.length))
        )
        dist[dist <= 0] <- 0


        if(length(random.length) == 1){
            random.length.p <-random.length
            random.length.n <-random.length
        }else{
            random.length.p <-random.length[p.index]
            random.length.n <-random.length[n.index]
        }

        # positive strand
        if(length(p.index) != 0){
            dist.p <- dist[p.index]
            regions.p <- regions.A[p.index]
            r.start.p <- min(start(regions.p))
            r.names.p <- names(regions.p)

            disp.p <- data.frame(start = r.start.p,
                                 distance = dist.p,
                                 names = r.names.p)
            R.p <- calculateShift(regions = regions.p,
                                  disp = disp.p,
                                  direction = 'right',
                                  strand = '+')

            disp.B.p <- data.frame(start = end(R.p),
                                   distance = random.length.p -1,
                                   names = r.names.p)

            B.p <- calculateShift(regions = regions.p,
                                  disp = disp.B.p,
                                  direction = 'right',
                                  strand = '+')

            B.p.exons <- extractRegions(regions_A = regions.p, B.p, strand = '+')
        }else{
            B.p.exons <- c()
        }

        # negative strand
        if(length(n.index) != 0){
            dist.n <- dist[n.index]
            regions.n <- regions.A[n.index]
            r.start.n <- max(end(regions.n))
            r.names.n <- names(regions.n)

            disp.n <- data.frame(start = r.start.n,
                                 distance = dist.n,
                                 names = r.names.n,
                                 strand = '-')
            R.n <- calculateShift(regions = regions.n,
                                  disp = disp.n,
                                  direction = 'right',
                                  strand = '-')

            disp.B.n <- data.frame(start = start(R.n),
                                   distance = random.length.n -1,
                                   names = r.names.n,
                                   strand = '-')
            B.n <- calculateShift(regions = regions.n,
                                  disp = disp.B.n,
                                  direction = 'right',
                                  strand = '-')

            B.n.exons <- extractRegions(regions_A = regions.n, B.n, strand = '-')
            B.n.exons$group <- B.n.exons$group + max(B.p.exons$group)
        }else{
            B.n.exons <-c()
        }
        randomResults <- c(B.p.exons, B.n.exons)
        randomResults <- GRanges2GRangesList(randomResults)
        return(randomResults)
    }
