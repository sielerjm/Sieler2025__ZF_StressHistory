# Simplified Neutral Community Model Function
# Based on Sloan et al. 2006 neutral model
# This version focuses on NLS fitting and avoids problematic MLE steps

sncm.fit.simple <- function(spp, pool = NULL, stats = TRUE, taxon = NULL) {
    require(minpack.lm)
    require(Hmisc)

    options(warn = -1)

    # Calculate the number of individuals per community
    N <- mean(apply(spp, 1, sum))

    # Calculate the average relative abundance of each taxa across communities
    if (is.null(pool)) {
        p.m <- apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m / N
    } else {
        p.m <- apply(pool, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m / N
    }

    # Calculate the occurrence frequency of each taxa across communities
    spp.bi <- 1 * (spp > 0)
    freq <- apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]

    # Combine
    C <- merge(p, freq, by = 0)
    C <- C[order(C[, 2]), ]
    C <- as.data.frame(C)
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ] # Removes rows with any zero
    p <- C.0[, 2]
    freq <- C.0[, 3]
    names(p) <- C.0[, 1]
    names(freq) <- C.0[, 1]

    # Calculate the limit of detection
    d <- 1 / N

    ## Fit model parameter m using Non-linear least squares (NLS)
    m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
    
    # Get migration rate
    m.coef <- coef(m.fit)
    
    # Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- pbeta(d, N * m.coef * p, N * m.coef * (1 - p), lower.tail = FALSE)
    Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq - freq.pred)^2) / (length(freq) - 1))

    # Calculate confidence intervals for predictions
    pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)

    ## Results
    if (stats == TRUE) {
        # Create simplified stats dataframe
        fitstats <- data.frame(
            m = m.coef,
            Rsqr = Rsqr,
            RMSE = RMSE,
            N = N,
            Samples = nrow(spp),
            Richness = length(p),
            Detect = d
        )
        return(fitstats)
    } else {
        # Return predictions
        A <- cbind(p, freq, freq.pred, pred.ci[, 2:3])
        A <- as.data.frame(A)
        colnames(A) <- c("p", "freq", "freq.pred", "pred.lwr", "pred.upr")
        
        if (is.null(taxon)) {
            B <- A[order(A[, 1]), ]
        } else {
            B <- merge(A, taxon, by = 0, all = TRUE)
            row.names(B) <- B[, 1]
            B <- B[, -1]
            B <- B[order(B[, 1]), ]
        }
        return(B)
    }
} 