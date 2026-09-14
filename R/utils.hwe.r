################################################################################
############## utils functions for gl.report.hwe ###############################
################################################################################

####################################### GenerateSamples #######################
#' Enumerates all genotype compositions (AA, AB, BB) of a sample of size n
#' @param n Sample size (number of individuals)
#' @return A matrix with columns AA, AB, BB; one row per composition
#' @noRd
GenerateSamples <- function(n = 5) {
    # generates all possible samples of size n.
    Res <- NULL
    for (i in 0:n) {
        AA <- i
        for (j in 0:(n - i)) {
            AB <- j
            BB <- (n - (AA + AB))
            sam <- c(AA, AB, BB)
            Res <- rbind(Res, sam)
        }
    }
    rownames(Res) <- 1:nrow(Res)
    colnames(Res) <- c("AA", "AB", "BB")
    return(Res)
}

############################################## CritSam #######################
#' Critical genotype samples for the HWE exact test (HardyWeinberg::HWExact)
#' @param n Sample size; Dpos TRUE for D > 0 (heterozygote excess) side;
#' alphalimit significance limit; pvaluetype passed to HWExact
#' @return list(Xn = critical genotype frequencies, Ds, fA)
#' @noRd
CritSam <- function(n, Dpos, alphalimit, pvaluetype) {
    X <- GenerateSamples(n)
    Res <- NULL
    Ds <- NULL
    pval <- NULL
    fA <- NULL
    for (i in 1:nrow(X)) {
        fA <- c(fA, (2 * X[i, 1] + X[i, 2]) / (2 * n))
        Ds <- suppressWarnings(
            c(Ds, HardyWeinberg::HWChisq(X[i,], verbose = FALSE)$D))
        pval <- suppressWarnings(
            c(
                pval,
                HardyWeinberg::HWExact(
                    X[i,],
                    alternative = "two.sided",
                    pvaluetype = pvaluetype,
                    verbose = FALSE
                )$pval
            ))
    }
    
    Y <- data.frame(X[, 1], X[, 2], X[, 3], fA, Ds, pval)
    colnames(Y) <- c("AA", "AB", "BB", "fA", "Ds", "pval")
    if (Dpos)
        Y <- Y[Y$Ds > 0,]
    else
        Y <- Y[Y$Ds < 0,]
    Y <- Y[Y$pval < alphalimit,]
    fre <- unique(fA)
    for (i in 1:length(fre)) {
        Ys <- Y[Y$fA == fre[i],]
        if (nrow(Ys) > 0) {
            indi <- which.max(Ys$pval)
            Ys <- Ys[indi,]
            Res <- rbind(Res, c(Ys$AA, Ys$AB, Ys$BB))
        }
    }
    Xn <- Res / n
    return(list(Xn = Xn, Ds = Ds, fA = fA))
}

############################################ CritSam_Chi #######################
#' Critical genotype samples for the HWE chi-square test
#' (HardyWeinberg::HWChisq); cc is the continuity correction
#' @return list(Xn = critical genotype frequencies, Ds, fA)
#' @noRd
CritSam_Chi <- function(n, Dpos, alphalimit, cc) {
    X <- GenerateSamples(n)
    Res <- NULL
    Ds <- NULL
    pval <- NULL
    fA <- NULL
    for (i in 1:nrow(X)) {
        fA <- c(fA, (2 * X[i, 1] + X[i, 2]) / (2 * n))
        Ds <- suppressWarnings(
            c(Ds,
              HardyWeinberg::HWChisq(X[i,], cc = cc, verbose = FALSE)$D))
        pval <- suppressWarnings(
            c(pval,
              HardyWeinberg::HWChisq(X[i,], cc = cc, verbose = FALSE)$pval))
    }
    
    Y <- data.frame(X[, 1], X[, 2], X[, 3], fA, Ds, pval)
    colnames(Y) <- c("AA", "AB", "BB", "fA", "Ds", "pval")
    if (Dpos)
        Y <- Y[Y$Ds > 0,]
    else
        Y <- Y[Y$Ds < 0,]
    Y <- Y[Y$pval < alphalimit,]
    fre <- unique(fA)
    for (i in 1:length(fre)) {
        Ys <- Y[Y$fA == fre[i],]
        if (nrow(Ys) > 0) {
            indi <- which.max(Ys$pval)
            Ys <- Ys[indi,]
            Res <- rbind(Res, c(Ys$AA, Ys$AB, Ys$BB))
        }
    }
    Xn <- Res / n
    return(list(Xn = Xn, Ds = Ds, fA = fA))
}
