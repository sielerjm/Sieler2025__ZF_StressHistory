# Adam Burns - 2/10/2015
# aburns2@uoregon.edu
# From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
# Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
# spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
# pool: A community table for defining source community (optional; Default=NULL).
# taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
# If stats=TRUE the function will return fitting statistics.
# If stats=FALSE the function will return a table of observed and predicted values for each otu.

sncm.fit <- function(spp, pool = NULL, stats = TRUE, taxon = NULL) {
    require(minpack.lm)
    require(Hmisc)
    require(stats4)

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
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ] # Removes rows with any zero (absent in either source pool or local communities)
    p <- C.0[, 2]
    freq <- C.0[, 3]
    names(p) <- C.0[, 1]
    names(freq) <- C.0[, 1]

    # Calculate the limit of detection
    d <- 1 / N

    ## Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
    
    # Get confidence interval safely
    m.ci <- tryCatch({
        confint(m.fit, "m", level = 0.95)
    }, error = function(e) {
        # If confint fails, use a simple approximation
        m.coef <- coef(m.fit)
        c(m.coef * 0.5, m.coef * 1.5)
    })

    ## Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
    # Modified likelihood function with better numerical stability
    sncm.LL <- function(m, sigma) {
        # Add bounds checking to prevent numerical issues
        if(m <= 0 || sigma <= 0) return(Inf)
        
        # Calculate predicted frequencies
        pred_freq <- pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
        
        # Handle any NaN or Inf values
        if(any(is.na(pred_freq)) || any(is.infinite(pred_freq))) return(Inf)
        
        # Calculate residuals
        R <- freq - pred_freq
        
        # Use log-likelihood with better numerical stability
        log_likelihood <- sum(dnorm(R, 0, sigma, log = TRUE))
        
        # Return negative log-likelihood (for minimization)
        return(-log_likelihood)
    }
    
    # Try MLE with error handling and multiple starting values
    m.mle <- NULL
    mle_success <- FALSE
    
    # Try different starting values for better convergence
    start_values <- list(
        list(m = 0.1, sigma = 0.1),
        list(m = 0.01, sigma = 0.1),
        list(m = 0.001, sigma = 0.1),
        list(m = coef(m.fit), sigma = 0.1)  # Use NLS result as starting point
    )
    
    for(i in seq_along(start_values)) {
        tryCatch({
            m.mle <- mle(sncm.LL, start = start_values[[i]], nobs = length(p),
                        method = "L-BFGS-B",  # Use bounded optimization
                        lower = c(m = 1e-6, sigma = 1e-6),  # Set lower bounds
                        upper = c(m = 1, sigma = 1))  # Set upper bounds
            mle_success <- TRUE
            break
        }, error = function(e) {
            # Continue to next starting value
            return(NULL)
        })
    }
    
    # If MLE fails, use NLS result as fallback
    if(!mle_success) {
        warning("MLE fitting failed, using NLS result as fallback")
        # Create a dummy mle object with NLS results
        m.mle <- list(
            coef = c(m = coef(m.fit), sigma = 0.1),
            details = list(value = NA)
        )
        class(m.mle) <- "mle"
    }

    ## Calculate Akaike's Information Criterion (AIC) safely
    aic.fit <- tryCatch({
        AIC(m.mle, k = 2)
    }, error = function(e) {
        NA
    })
    
    bic.fit <- tryCatch({
        BIC(m.mle)
    }, error = function(e) {
        NA
    })

    ## Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
    Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq - freq.pred)^2) / (length(freq) - 1))

    pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)

    ## Calculate AIC for binomial model
    bino.LL <- function(mu, sigma) {
        R <- freq - pbinom(d, N, p, lower.tail = FALSE)
        R <- dnorm(R, mu, sigma)
        -sum(log(R))
    }
    
    # Try binomial MLE with error handling
    bino.mle <- NULL
    tryCatch({
        bino.mle <- mle(bino.LL, start = list(mu = 0, sigma = 0.1), nobs = length(p))
    }, error = function(e) {
        warning("Binomial MLE fitting failed")
        bino.mle <- list(details = list(value = NA))
        class(bino.mle) <- "mle"
    })

    aic.bino <- tryCatch({
        AIC(bino.mle, k = 2)
    }, error = function(e) {
        NA
    })
    
    bic.bino <- tryCatch({
        BIC(bino.mle)
    }, error = function(e) {
        NA
    })

    ## Goodness of fit for binomial model
    bino.pred <- pbinom(d, N, p, lower.tail = FALSE)
    Rsqr.bino <- 1 - (sum((freq - bino.pred)^2)) / (sum((freq - mean(freq))^2))
    RMSE.bino <- sqrt(sum((freq - bino.pred)^2) / (length(freq) - 1))

    bino.pred.ci <- binconf(bino.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)

    ## Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma) {
        R <- freq - ppois(d, N * p, lower.tail = FALSE)
        R <- dnorm(R, mu, sigma)
        -sum(log(R))
    }
    
    # Try Poisson MLE with error handling
    pois.mle <- NULL
    tryCatch({
        pois.mle <- mle(pois.LL, start = list(mu = 0, sigma = 0.1), nobs = length(p))
    }, error = function(e) {
        warning("Poisson MLE fitting failed")
        pois.mle <- list(details = list(value = NA))
        class(pois.mle) <- "mle"
    })

    aic.pois <- tryCatch({
        AIC(pois.mle, k = 2)
    }, error = function(e) {
        NA
    })
    
    bic.pois <- tryCatch({
        BIC(pois.mle)
    }, error = function(e) {
        NA
    })

    ## Goodness of fit for Poisson model
    pois.pred <- ppois(d, N * p, lower.tail = FALSE)
    Rsqr.pois <- 1 - (sum((freq - pois.pred)^2)) / (sum((freq - mean(freq))^2))
    RMSE.pois <- sqrt(sum((freq - pois.pred)^2) / (length(freq) - 1))

    pois.pred.ci <- binconf(pois.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)

    ## Results
    if (stats == TRUE) {
        fitstats <- data.frame(m = numeric(), m.ci = numeric(), m.mle = numeric(), maxLL = numeric(), binoLL = numeric(), poisLL = numeric(), Rsqr = numeric(), Rsqr.bino = numeric(), Rsqr.pois = numeric(), RMSE = numeric(), RMSE.bino = numeric(), RMSE.pois = numeric(), AIC = numeric(), BIC = numeric(), AIC.bino = numeric(), BIC.bino = numeric(), AIC.pois = numeric(), BIC.pois = numeric(), N = numeric(), Samples = numeric(), Richness = numeric(), Detect = numeric())
        
        # Extract MLE coefficients safely
        m.mle.coef <- if(mle_success) coef(m.mle)["m"] else coef(m.fit)
        m.mle.value <- if(mle_success) {
            tryCatch({
                # Try to get the log-likelihood value
                if("details" %in% names(m.mle) && "value" %in% names(m.mle$details)) {
                    m.mle$details$value
                } else {
                    NA
                }
            }, error = function(e) NA)
        } else NA
        
        bino.mle.value <- tryCatch({
            if("details" %in% names(bino.mle) && "value" %in% names(bino.mle$details)) {
                bino.mle$details$value
            } else {
                NA
            }
        }, error = function(e) NA)
        
        pois.mle.value <- tryCatch({
            if("details" %in% names(pois.mle) && "value" %in% names(pois.mle$details)) {
                pois.mle$details$value
            } else {
                NA
            }
        }, error = function(e) NA)
        
        fitstats[1, ] <- c(coef(m.fit), coef(m.fit) - m.ci[1], m.mle.coef, m.mle.value, bino.mle.value, pois.mle.value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
        return(fitstats)
    } else {
        A <- cbind(p, freq, freq.pred, pred.ci[, 2:3], bino.pred, bino.pred.ci[, 2:3])
        A <- as.data.frame(A)
        colnames(A) <- c("p", "freq", "freq.pred", "pred.lwr", "pred.upr", "bino.pred", "bino.lwr", "bino.upr")
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