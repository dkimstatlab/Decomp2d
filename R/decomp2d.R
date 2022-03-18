# version 0.6.0

empperiod <- function(z) {
    rc <- dim(z)
    rowperiod <- colperiod <- NULL
    for (i in 1:rc[2])
        rowperiod <- c(rowperiod, diff(localextrema(z[,i])$maxindex[,1]))

    for (i in 1:rc[1])
        colperiod <- c(colperiod, diff(localextrema(z[i,])$maxindex[,1]))

    list(rowperiod=rowperiod, colperiod=colperiod)
}

decomp2d <- function(z, method="wavelet",
                     dct.frequency=NULL,
                     emd.sm=FALSE, emd.spar=NULL, emd.tol = 0.1^2, emd.maxiter = 20,
                     ept.tau=NULL, ept.process = c("average", "average"), ept.tol=0.1^2, ept.maxiter = 20,
                     pca.freqcomp=NULL,
                     ssa.L=NULL, ssa.freqcomp=NULL,
                     wavelet.highlevel=NULL) {
    if (method == "ept") {
        if (is.null(ept.tau))
            stop ("It is required to specify a size parameter for two-dimensional ensemble patch transform
                  when ept is used." )
        eptz <- eptdecomp2d(z=z, type="rectangle", tau=ept.tau, process = ept.process, tol=ept.tol, maxiter = ept.maxiter, check = FALSE)
        residue <- eptz$residue

    }
    if (method == "emd") {
        sm <- ifelse(emd.sm, "locfit", "none")
        if (emd.sm & is.null(emd.spar))
            stop ("It is required to specify user-supplied smoothing parameter of local polynomial smoothing
                  when smoothing emd is used." )
        emdz <- extractimf2d(z, tol = emd.tol, max.sift = emd.maxiter, boundary = "reflexive", sm=sm, spar=emd.spar, check = FALSE)
        residue <- emdz$residue
    }
    if (method == "dct") {
        if (is.null(dct.frequency))
            stop ("It is required to specify user-supplied threshold of frequencies
                  when dct is used." )
        dctz <- DCT2D(z, returnmat=TRUE)

        # if (!is.null(dct.quant)) {
        #     thresh <- quantile(abs(c(dctz)), dct.quant, type=3)
        #     dct.frequency <- which(thresh == abs(dctz), arr.ind=TRUE)
        #     lowcomp <- matrix(0, nrow=dim(dctz)[1], ncol=dim(dctz)[1])
        #     lowcomp[1:dct.frequency[1], 1:dct.frequency[2]] <- dctz[1:dct.frequency[1], 1:dct.frequency[2]]
        #     residue <- IDCT2D(lowcomp, returnmat=TRUE)
        # }
        # if (!is.null(dct.frequency)) {
            if (length(dct.frequency) == 1) dct.frequency <- rep(dct.frequency, 2)
            lowcomp <- matrix(0, nrow=dim(dctz)[1], ncol=dim(dctz)[1])
            lowcomp[1:dct.frequency[1], 1:dct.frequency[2]] <- dctz[1:dct.frequency[1], 1:dct.frequency[2]]
            residue <- IDCT2D(lowcomp, returnmat=TRUE)
        # }
    }
    if (method == "ssa") {
        if (is.null(ssa.L) | is.null(ssa.freqcomp))
            stop ("It is required to specify user-supplied window length and numeric vectors of frequency components
                  when ssa is used." )
        ssaz <- ssa(z, kind="2d-ssa", L=ssa.L)
        residue <- z - reconstruct(ssaz, groups=list(fc=ssa.freqcomp))$fc
    }
    if (method == "pca") {
        if (is.null(pca.freqcomp))
            stop ("It is required to specify user-supplied numeric vectors of frequency components
                  when pca is used.")
        pcaz <- svd(t(z - mean(z)) %*% (z - mean(z)))
        pca.freqcomp <- sort(pca.freqcomp)
        residue <- 0
        for (i in 1:length(pca.freqcomp))
            residue <- residue + z %*% outer(pcaz$u[, pca.freqcomp[i]], pcaz$u[, pca.freqcomp[i]])
        residue <- z - residue
    }
    if (method == "wavelet") {
        if (is.null(wavelet.highlevel))
          stop ("It is required to specify a particular resolution level of high-frequency component
                      when wavelet decomposotion is used." )
        imwdz <- imwd(z, bc="symmetric")
        imwdthresh <- nullevels(imwdz , levelstonull=wavelet.highlevel:(nlevelsWT(imwdz)-1))
        imwdrecon <- imwr(imwdthresh)
        residue <- imwdrecon
    }
    list(fc=z-residue, residue=residue)
}


