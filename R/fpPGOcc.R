fpPGOcc <- function(occ.formula, tp.det.formula, fp.det.formula, 
                    data, inits, priors, tuning, n.batch, batch.length, 
                    accept.rate = 0.43, n.omp.threads = 1, verbose = TRUE, 
                    n.report = 100, n.burn = round(.10 * n.batch * batch.length),
                    n.thin = 1, n.chains = 1, ...) {
  
  ptm <- proc.time()
  
  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n");
    cat("\tPreparing to run the model\n");
    cat("----------------------------------------\n");
  }
  
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  # Call ----------------------------------------------------------------
  # Returns a call in which all of the specified arguments are
  # specified by their full names.
  cl <- match.call()
  
  # Some initial checks -------------------------------------------------
  if (missing(data)) {
    stop("data must be specified")
  }
  if (!is.list(data)) {
    stop("data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(occ.formula)) {
    stop("occ.formula must be specified")
  }
  if (missing(tp.det.formula)) {
    stop("tp.det.formula must be specified")
  }
  if (missing(fp.det.formula)) {
    stop("fp.det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("detection data y must be specified in data")
  }
  y <- as.matrix(data$y)
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, dim(y)[1], 1)
    } else {
      stop("occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!is.matrix(data$occ.covs) & !is.data.frame(data$occ.covs)) {
    stop("occ.covs must be a matrix or data frame")
  }
  if (dim(y)[[1]] != dim(data$occ.covs)[[1]]){
    stop("number of sites in detection data 'y' (J = ", dim(y)[[1]],
         ") does not match number of sites in site-specific covariates 'occ.covs' (J = ", 
         dim(data$occ.covs)[[1]], ").")
  }
  
  # Check validation data
  if (!'n.val' %in% names(data)) {
    stop("human validation data n.val must be specified in data")
  }
  if (!'d' %in% names(data)) {
    stop("true positive human validation data d must be specified in data")
  }
  n.val <- as.matrix(data$n.val)
  d <- as.matrix(data$d)
  if (sum(dim(y) != dim(n.val)) > 0) {
    stop("dimension mismatch: number of sites and/or replicates does not match between y and n.val")
  }
  if (sum(dim(y) != dim(d)) > 0) {
    stop("dimension mismatch: number of sites and/or replicates does not match between y and d")
  }
  
  if (sum(is.na(data$occ.covs)) > 0) {
    stop("missing covariate values in data$occ.covs. Remove these sites from all data or impute non-missing values.")
  }
  if (!'tp.det.covs' %in% names(data)) {
    if (tp.det.formula == ~ 1) {
      if (verbose) {
        message("true positive detection covariates (tp.det.covs) not specified in data.\nAssuming intercept only true positive detection model.\n")
}
      data$tp.det.covs <- list(int = rep(1, dim(y)[1]))
    } else {
      stop("tp.det.covs must be specified in data for a true positive detection model with covariates")
    }
  }
  if (!is.list(data$tp.det.covs)) {
    stop("tp.det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (!'fp.det.covs' %in% names(data)) {
    if (fp.det.formula == ~ 1) {
      if (verbose) {
        message("false positive detection covariates (fp.det.covs) not specified in data.\nAssuming intercept only false positive detection model.\n")
}
      data$fp.det.covs <- list(int = rep(1, dim(y)[1]))
    } else {
      stop("fp.det.covs must be specified in data for a true positive detection model with covariates")
    }
  }
  if (!is.list(data$fp.det.covs)) {
    stop("fp.det.covs must be a list of matrices, data frames, and/or vectors")
  }
  
  if (missing(n.batch)) {
    stop("must specify number of MCMC batches")
  }
  if (missing(batch.length)) {
    stop("must specify length of each MCMC batch")
  }
  n.samples <- n.batch * batch.length
  if (n.burn > n.samples) {
    stop("n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("n.thin must be less than n.samples")
  }
  # Check if n.burn, n.thin, and n.samples result in an integer and error if otherwise.
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number. Please respecify the MCMC criteria such that the number of posterior samples saved is a whole number.")
  }
  
  # First subset all detection covariates to only use those that are included in the analysis. 
  data$tp.det.covs <- data$tp.det.covs[names(data$tp.det.covs) %in% all.vars(tp.det.formula)]
  data$fp.det.covs <- data$fp.det.covs[names(data$fp.det.covs) %in% all.vars(fp.det.formula)]
  # Null model support
  if (length(data$tp.det.covs) == 0) {
    data$tp.det.covs <- list(int = rep(1, dim(y)[1]))
  }
  if (length(data$fp.det.covs) == 0) {
    data$fp.det.covs <- list(int = rep(1, dim(y)[1]))
  }
  # Make both covariates a data frame. Unlist is necessary for when factors
  # are supplied. 
  data$tp.det.covs <- data.frame(lapply(data$tp.det.covs, function(a) unlist(c(a))))
  data$fp.det.covs <- data.frame(lapply(data$fp.det.covs, function(a) unlist(c(a))))
  # Check if all detection covariates are at site level, and simplify the data
  # if necessary. This is highly unlikely in practice
  y.big <- y
  n.val.big <- n.val
  d.big <- d
  
  
  data$occ.covs <- as.data.frame(data$occ.covs)
  
  # Check whether random effects are sent in as numeric, and
  # return error if they are not. 
  # Occurrence ----------------------
  if (!is.null(findbars(occ.formula))) {
    occ.re.names <- sapply(findbars(occ.formula), all.vars)
    for (i in 1:length(occ.re.names)) {
      if (is(data$occ.covs[, occ.re.names[i]], 'factor')) {
        stop(paste("random effect variable ", occ.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$occ.covs[, occ.re.names[i]], 'character')) {
        stop(paste("random effect variable ", occ.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # TP Detection ----------------------
  if (!is.null(findbars(tp.det.formula))) {
    tp.det.re.names <- sapply(findbars(tp.det.formula), all.vars)
    for (i in 1:length(tp.det.re.names)) {
      if (is(data$tp.det.covs[, tp.det.re.names[i]], 'factor')) {
        stop(paste("random effect variable ", tp.det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$tp.det.covs[, tp.det.re.names[i]], 'character')) {
        stop(paste("random effect variable ", tp.det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  # FP Detection ----------------------
  if (!is.null(findbars(fp.det.formula))) {
    fp.det.re.names <- sapply(findbars(fp.det.formula), all.vars)
    for (i in 1:length(fp.det.re.names)) {
      if (is(data$fp.det.covs[, fp.det.re.names[i]], 'factor')) {
        stop(paste("random effect variable ", fp.det.re.names[i], " specified as a factor. Random effect variables must be specified as numeric.", sep = ''))
      } 
      if (is(data$fp.det.covs[, fp.det.re.names[i]], 'character')) {
        stop(paste("random effect variable ", fp.det.re.names[i], " specified as character. Random effect variables must be specified as numeric.", sep = ''))
      }
    }
  }
  
  # Checking missing values ---------------------------------------------
  # y -------------------------------
  y.na.test <- apply(y.big, 1, function(a) sum(!is.na(a)))
  if (sum(y.na.test == 0) > 0) {
    stop("some sites in y have all missing detection histories. Remove these sites from all objects in the 'data' argument, then use 'predict' to obtain predictions at these locations if desired.")
  }
  # occ.covs ------------------------
  if (sum(is.na(data$occ.covs)) != 0) {
    stop("missing values in occ.covs. Please remove these sites from all objects in data or somehow replace the NA values with non-missing values (e.g., mean imputation).") 
  }
  # tp.det.covs ------------------------
  for (i in 1:ncol(data$tp.det.covs)) {
    if (sum(is.na(data$tp.det.covs[, i])) > sum(is.na(y.big))) {
      stop("some elements in tp.det.covs have missing values where there is an observed data value in y. Please either replace the NA values in tp.det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  # Misalignment between y and det.covs
  y.missing <- which(is.na(data$y))
  tp.det.covs.missing <- lapply(data$tp.det.covs, function(a) which(is.na(a)))
  for (i in 1:length(tp.det.covs.missing)) {
    tmp.indx <- !(y.missing %in% tp.det.covs.missing[[i]])
    if (sum(tmp.indx) > 0) {
      if (i == 1 & verbose) {
        message("There are missing values in data$y with corresponding non-missing values in data$tp.det.covs.\nRemoving these site/replicate combinations for fitting the model.")
      }
      data$tp.det.covs[y.missing, i] <- NA
    }
  }
  # fp.det.covs ------------------------
  for (i in 1:ncol(data$fp.det.covs)) {
    if (sum(is.na(data$fp.det.covs[, i])) > sum(is.na(y.big))) {
      stop("some elements in fp.det.covs have missing values where there is an observed data value in y. Please either replace the NA values in fp.det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  # Misalignment between y and det.covs
  y.missing <- which(is.na(data$y))
  fp.det.covs.missing <- lapply(data$fp.det.covs, function(a) which(is.na(a)))
  for (i in 1:length(fp.det.covs.missing)) {
    tmp.indx <- !(y.missing %in% fp.det.covs.missing[[i]])
    if (sum(tmp.indx) > 0) {
      if (i == 1 & verbose) {
        message("There are missing values in data$y with corresponding non-missing values in data$fp.det.covs.\nRemoving these site/replicate combinations for fitting the model.")
      }
      data$fp.det.covs[y.missing, i] <- NA
    }
  }
  # fp.det.covs ------------------------
  for (i in 1:ncol(data$fp.det.covs)) {
    if (sum(is.na(data$fp.det.covs[, i])) > sum(is.na(y.big))) {
      stop("some elements in fp.det.covs have missing values where there is an observed data value in y. Please either replace the NA values in fp.det.covs with non-missing values (e.g., mean imputation) or set the corresponding values in y to NA where the covariate is missing.") 
    }
  }
  # Misalignment between y and det.covs
  y.missing <- which(is.na(data$y))
  fp.det.covs.missing <- lapply(data$fp.det.covs, function(a) which(is.na(a)))
  for (i in 1:length(fp.det.covs.missing)) {
    tmp.indx <- !(y.missing %in% fp.det.covs.missing[[i]])
    if (sum(tmp.indx) > 0) {
      if (i == 1 & verbose) {
        message("There are missing values in data$y with corresponding non-missing values in data$fp.det.covs.\nRemoving these site/replicate combinations for fitting the model.")
      }
      data$fp.det.covs[y.missing, i] <- NA
    }
  }
  # TODO: 
  # Check missingness patterns in the validated data
  # Cells that are missing should be identical between d, n.val, and y. 
  n.val.missing <- which(is.na(data$n.val)) 
  d.missing <- which(is.na(data$d))
  if (!all.equal(n.val.missing, y.missing)) {
    stop("mismatch between missing values in n.val and y. These must have the same number and location of missing values")
  }
  if (!all.equal(n.val.missing, d.missing)) {
    stop("mismatch between missing values in n.val and d. These must have the same number and location of missing values")
  }
  
  # Formula -------------------------------------------------------------
  # Occupancy -----------------------
  if (is(occ.formula, 'formula')) {
    # NOTE: this will reorder random effect names by the number of unique 
    #       levels they have, not necessarily the order in which they are specified
    #       in the model formula. 
    tmp <- parseFormula(occ.formula, data$occ.covs)
    X <- as.matrix(tmp[[1]])
    X.re <- as.matrix(tmp[[4]])
    x.re.names <- colnames(X.re)
    x.names <- tmp[[2]]
  } else {
    stop("occ.formula is misspecified")
  }
  # Get RE level names
  re.level.names <- lapply(data$occ.covs[, x.re.names, drop = FALSE],
                           function (a) sort(unique(a)))
  
  # TP Detection ---------------------
  if (is(tp.det.formula, 'formula')) {
    tmp <- parseFormula(tp.det.formula, data$tp.det.covs)
    X.mu <- as.matrix(tmp[[1]])
    X.mu.re <- as.matrix(tmp[[4]])
    x.mu.re.names <- colnames(X.mu.re)
    x.mu.names <- tmp[[2]]
  } else {
    stop("tp.det.formula is misspecified")
  }
  mu.re.level.names <- lapply(data$tp.det.covs[, x.mu.re.names, drop = FALSE],
                             function (a) sort(unique(a)))
  # FP Detection ---------------------
  if (is(fp.det.formula, 'formula')) {
    tmp <- parseFormula(fp.det.formula, data$fp.det.covs)
    X.omega <- as.matrix(tmp[[1]])
    X.omega.re <- as.matrix(tmp[[4]])
    x.omega.re.names <- colnames(X.omega.re)
    x.omega.names <- tmp[[2]]
  } else {
    stop("fp.det.formula is misspecified")
  }
  omega.re.level.names <- lapply(data$fp.det.covs[, x.omega.re.names, drop = FALSE],
                                 function (a) sort(unique(a)))
  
  # Get basic info from inputs ------------------------------------------
  # Number of sites
  J <- nrow(y.big)
  # Number of occupancy parameters
  p.occ <- ncol(X)
  # Number of occupancy random effect parameters
  p.occ.re <- ncol(X.re)
  # Number of tp detection parameters
  p.tp.det <- ncol(X.mu)
  # Number of tp detection random effect parameters
  p.tp.det.re <- ncol(X.mu.re)
  # Number of fp detection parameters
  p.fp.det <- ncol(X.omega)
  # Number of fp detection random effect parameters
  p.fp.det.re <- ncol(X.omega.re)
   
  # Number of latent occupancy random effect values
  n.occ.re <- length(unlist(apply(X.re, 2, unique)))
  n.occ.re.long <- apply(X.re, 2, function(a) length(unique(a)))
  # Number of latent tp detection random effect values
  n.tp.det.re <- length(unlist(apply(X.mu.re, 2, unique)))
  n.tp.det.re.long <- apply(X.mu.re, 2, function(a) length(unique(a)))
  if (p.tp.det.re == 0) n.tp.det.re.long <- 0
  # Number of latent fp detection random effect values
  n.fp.det.re <- length(unlist(apply(X.omega.re, 2, unique)))
  n.fp.det.re.long <- apply(X.omega.re, 2, function(a) length(unique(a)))
  if (p.fp.det.re == 0) n.fp.det.re.long <- 0
  # Number of replicates at each site
  n.rep <- apply(y.big, 1, function(a) sum(!is.na(a)))
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- which(!is.na(y.big[j, ]))
  }
  # Max number of repeat visits
  K.max <- max(n.rep)
  # Because I like K better than n.rep
  K <- n.rep
  
  # Get indices to map z to y -------------------------------------------
  z.long.indx <- rep(1:J, dim(y.big)[2])
  z.long.indx <- z.long.indx[!is.na(c(y.big))]
  # Subtract 1 for indices in C
  z.long.indx <- z.long.indx - 1
  # y, n.val, d are stored in the following order: site, visit
  y <- c(y)
  names.long <- which(!is.na(y))
  n.val <- c(n.val)
  d <- c(d)
  # Remove missing observations when the covariate data are available but
  # there are missing detection-nondetection data. 
  if (nrow(X.mu) == length(y)) {
    X.mu <- X.mu[!is.na(y), , drop = FALSE]  
  }
  if (nrow(X.omega) == length(y)) {
    X.omega <- X.omega[!is.na(y), , drop = FALSE]  
  }
  y <- y[!is.na(y)]
  n.val <- n.val[!is.na(n.val)]
  d <- d[!is.na(d)]
  # Number of data points for the y vector
  n.obs <- nrow(X.mu)
  
  # Get random effect matrices all set ----------------------------------
  if (p.occ.re > 1) {
    for (j in 2:p.occ.re) {
      X.re[, j] <- X.re[, j] + max(X.re[, j - 1]) + 1
    }
  }
  if (p.tp.det.re > 1) {
    for (j in 2:p.tp.det.re) {
      X.mu.re[, j] <- X.mu.re[, j] + max(X.mu.re[, j - 1]) + 1
    }
  }
  if (p.fp.det.re > 1) {
    for (j in 2:p.fp.det.re) {
      X.omega.re[, j] <- X.omega.re[, j] + max(X.omega.re[, j - 1]) + 1
    }
  }

  # Priors --------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))
  # beta -----------------------
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p.occ & length(mu.beta) != 1) {
      if (p.occ == 1) {
        stop(paste("beta.normal[[1]] must be a vector of length ",
                   p.occ, " with elements corresponding to betas' mean", sep = ""))
      } else {
        stop(paste("beta.normal[[1]] must be a vector of length ",
                   p.occ, " or 1 with elements corresponding to betas' mean", sep = ""))
      }
    }
    if (length(sigma.beta) != p.occ & length(sigma.beta) != 1) {
      if (p.occ == 1) {
        stop(paste("beta.normal[[2]] must be a vector of length ",
                   p.occ, " with elements corresponding to betas' variance", sep = ""))
      } else {
        stop(paste("beta.normal[[2]] must be a vector of length ",
                   p.occ, " or 1 with elements corresponding to betas' variance", sep = ""))
      }
    }
    if (length(sigma.beta) != p.occ) {
      sigma.beta <- rep(sigma.beta, p.occ)
    }
    if (length(mu.beta) != p.occ) {
      mu.beta <- rep(mu.beta, p.occ)
    }
    Sigma.beta <- sigma.beta * diag(p.occ)
  } else {
    if (verbose) {
      message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.beta <- rep(0, p.occ)
    sigma.beta <- rep(2.72, p.occ)
    Sigma.beta <- diag(p.occ) * 2.72
  }
  # gamma -----------------------
  if ("gamma.normal" %in% names(priors)) {
    if (!is.list(priors$gamma.normal) | length(priors$gamma.normal) != 2) {
      stop("gamma.normal must be a list of length 2")
    }
    mu.gamma <- priors$gamma.normal[[1]]
    sigma.gamma <- priors$gamma.normal[[2]]
    if (length(mu.gamma) != p.tp.det & length(mu.gamma) != 1) {
      if (p.tp.det == 1) {
        stop(paste("gamma.normal[[1]] must be a vector of length ",
                   p.tp.det, " with elements corresponding to gammas' mean", sep = ""))
      } else {
        stop(paste("gamma.normal[[1]] must be a vector of length ",
                   p.tp.det, " or 1 with elements corresponding to gammas' mean", sep = ""))
      }
    }
    if (length(sigma.gamma) != p.tp.det & length(sigma.gamma) != 1) {
      if (p.tp.det == 1) {
        stop(paste("gamma.normal[[2]] must be a vector of length ",
                   p.tp.det, " with elements corresponding to gammas' variance", sep = ""))
      } else {
        stop(paste("gamma.normal[[2]] must be a vector of length ",
                   p.tp.det, " or 1 with elements corresponding to gammas' variance", sep = ""))
      }
    }
    if (length(sigma.gamma) != p.tp.det) {
      sigma.gamma <- rep(sigma.gamma, p.tp.det)
    }
    if (length(mu.gamma) != p.tp.det) {
      mu.gamma <- rep(mu.gamma, p.tp.det)
    }
    Sigma.gamma <- sigma.gamma * diag(p.tp.det)
  } else {
    if (verbose) {
      message("No prior specified for gamma.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.gamma <- rep(0, p.tp.det)
    sigma.gamma <- rep(2.72, p.tp.det)
    Sigma.gamma <- diag(p.tp.det) * 2.72
  }
  # sigma.sq.psi --------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.psi.ig) | length(priors$sigma.sq.psi.ig) != 2) {
        stop("sigma.sq.psi.ig must be a list of length 2")
      }
      sigma.sq.psi.a <- priors$sigma.sq.psi.ig[[1]]
      sigma.sq.psi.b <- priors$sigma.sq.psi.ig[[2]]
      if (length(sigma.sq.psi.a) != p.occ.re & length(sigma.sq.psi.a) != 1) {
        if (p.occ.re == 1) {
          stop(paste("sigma.sq.psi.ig[[1]] must be a vector of length ", 
                     p.occ.re, " with elements corresponding to sigma.sq.psis' shape", sep = ""))
        } else {
          stop(paste("sigma.sq.psi.ig[[1]] must be a vector of length ", 
                     p.occ.re, " or 1 with elements corresponding to sigma.sq.psis' shape", sep = ""))
        }
      }
      if (length(sigma.sq.psi.b) != p.occ.re & length(sigma.sq.psi.b) != 1) {
        if (p.occ.re == 1) {
          stop(paste("sigma.sq.psi.ig[[2]] must be a vector of length ", 
                     p.occ.re, " with elements corresponding to sigma.sq.psis' scale", sep = ""))
        } else {
          stop(paste("sigma.sq.psi.ig[[2]] must be a vector of length ", 
                     p.occ.re, " or 1with elements corresponding to sigma.sq.psis' scale", sep = ""))
        }
      }
      if (length(sigma.sq.psi.a) != p.occ.re) {
        sigma.sq.psi.a <- rep(sigma.sq.psi.a, p.occ.re)
      }
      if (length(sigma.sq.psi.b) != p.occ.re) {
        sigma.sq.psi.b <- rep(sigma.sq.psi.b, p.occ.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.psi.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.psi.a <- rep(0.1, p.occ.re)
      sigma.sq.psi.b <- rep(0.1, p.occ.re)
    }
  } else {
    sigma.sq.psi.a <- 0
    sigma.sq.psi.b <- 0
  }
  # delta -----------------------
  if ("delta.normal" %in% names(priors)) {
    if (!is.list(priors$delta.normal) | length(priors$delta.normal) != 2) {
      stop("delta.normal must be a list of length 2")
    }
    mu.delta <- priors$delta.normal[[1]]
    sigma.delta <- priors$delta.normal[[2]]
    if (length(mu.delta) != p.fp.det & length(mu.delta) != 1) {
      if (p.fp.det == 1) {
        stop(paste("delta.normal[[1]] must be a vector of length ",
                   p.fp.det, " with elements corresponding to deltas' mean", sep = ""))
      } else {
        stop(paste("delta.normal[[1]] must be a vector of length ",
                   p.fp.det, " or 1 with elements corresponding to deltas' mean", sep = ""))
      }
    }
    if (length(sigma.delta) != p.fp.det & length(sigma.delta) != 1) {
      if (p.fp.det == 1) {
        stop(paste("delta.normal[[2]] must be a vector of length ",
                   p.fp.det, " with elements corresponding to deltas' variance", sep = ""))
      } else {
        stop(paste("delta.normal[[2]] must be a vector of length ",
                   p.fp.det, " or 1 with elements corresponding to deltas' variance", sep = ""))
      }
    }
    if (length(sigma.delta) != p.fp.det) {
      sigma.delta <- rep(sigma.delta, p.fp.det)
    }
    if (length(mu.delta) != p.fp.det) {
      mu.delta <- rep(mu.delta, p.fp.det)
    }
    Sigma.delta <- sigma.delta * diag(p.fp.det)
  } else {
    if (verbose) {
      message("No prior specified for delta.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.delta <- rep(0, p.fp.det)
    sigma.delta <- rep(2.72, p.fp.det)
    Sigma.delta <- diag(p.fp.det) * 2.72
  }
  # sigma.sq.mu --------------------
  if (p.tp.det.re > 0) {
    if ("sigma.sq.mu.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.mu.ig) | length(priors$sigma.sq.mu.ig) != 2) {
        stop("sigma.sq.mu.ig must be a list of length 2")
      }
      sigma.sq.mu.a <- priors$sigma.sq.mu.ig[[1]]
      sigma.sq.mu.b <- priors$sigma.sq.mu.ig[[2]]
      if (length(sigma.sq.mu.a) != p.tp.det.re & length(sigma.sq.mu.a) != 1) {
        if (p.tp.det.re == 1) {
          stop(paste("sigma.sq.mu.ig[[1]] must be a vector of length ", 
                     p.tp.det.re, " with elements corresponding to sigma.sq.mu's shape", sep = ""))
        } else {
          stop(paste("sigma.sq.mu.ig[[1]] must be a vector of length ", 
                     p.tp.det.re, " or 1 with elements corresponding to sigma.sq.mu's shape", sep = ""))
        }
      }
      if (length(sigma.sq.mu.b) != p.tp.det.re & length(sigma.sq.mu.b) != 1) {
        if (p.tp.det.re == 1) {
          stop(paste("sigma.sq.mu.ig[[2]] must be a vector of length ", 
                     p.tp.det.re, " with elements corresponding to sigma.sq.mu's scale", sep = ""))
        } else {
          stop(paste("sigma.sq.mu.ig[[2]] must be a vector of length ", 
                     p.tp.det.re, " or 1 with elements corresponding to sigma.sq.mu's scale", sep = ""))
        }
      }
      if (length(sigma.sq.mu.a) != p.tp.det.re) {
        sigma.sq.mu.a <- rep(sigma.sq.mu.a, p.tp.det.re)
      }
      if (length(sigma.sq.mu.b) != p.tp.det.re) {
        sigma.sq.mu.b <- rep(sigma.sq.mu.b, p.tp.det.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.mu.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.mu.a <- rep(0.1, p.tp.det.re)
      sigma.sq.mu.b <- rep(0.1, p.tp.det.re)
    }
  } else {
    sigma.sq.mu.a <- 0
    sigma.sq.mu.b <- 0
  }
  # sigma.sq.omega --------------------
  if (p.fp.det.re > 0) {
    if ("sigma.sq.omega.ig" %in% names(priors)) {
      if (!is.list(priors$sigma.sq.omega.ig) | length(priors$sigma.sq.omega.ig) != 2) {
        stop("sigma.sq.omega.ig must be a list of length 2")
      }
      sigma.sq.omega.a <- priors$sigma.sq.omega.ig[[1]]
      sigma.sq.omega.b <- priors$sigma.sq.omega.ig[[2]]
      if (length(sigma.sq.omega.a) != p.fp.det.re & length(sigma.sq.omega.a) != 1) {
        if (p.fp.det.re == 1) {
          stop(paste("sigma.sq.omega.ig[[1]] must be a vector of length ", 
                     p.fp.det.re, " with elements corresponding to sigma.sq.omega's shape", sep = ""))
        } else {
          stop(paste("sigma.sq.omega.ig[[1]] must be a vector of length ", 
                     p.fp.det.re, " or 1 with elements corresponding to sigma.sq.omega's shape", sep = ""))
        }
      }
      if (length(sigma.sq.omega.b) != p.fp.det.re & length(sigma.sq.omega.b) != 1) {
        if (p.fp.det.re == 1) {
          stop(paste("sigma.sq.omega.ig[[2]] must be a vector of length ", 
                     p.fp.det.re, " with elements corresponding to sigma.sq.omega's scale", sep = ""))
        } else {
          stop(paste("sigma.sq.omega.ig[[2]] must be a vector of length ", 
                     p.fp.det.re, " or 1 with elements corresponding to sigma.sq.omega's scale", sep = ""))
        }
      }
      if (length(sigma.sq.omega.a) != p.fp.det.re) {
        sigma.sq.omega.a <- rep(sigma.sq.omega.a, p.fp.det.re)
      }
      if (length(sigma.sq.omega.b) != p.fp.det.re) {
        sigma.sq.omega.b <- rep(sigma.sq.omega.b, p.fp.det.re)
      }
    }   else {
      if (verbose) {	    
        message("No prior specified for sigma.sq.omega.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
      }
      sigma.sq.omega.a <- rep(0.1, p.fp.det.re)
      sigma.sq.omega.b <- rep(0.1, p.fp.det.re)
    }
  } else {
    sigma.sq.omega.a <- 0
    sigma.sq.omega.b <- 0
  }

  # Starting values -----------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z -------------------------------
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (!is.vector(z.inits)) {
      stop(paste("initial values for z must be a vector of length ",
                 J, sep = ""))
    }
    if (length(z.inits) != J) {
      stop(paste("initial values for z must be a vector of length ",
                 J, sep = ""))
    }
  } else {
    # Default is to set initial values based on wherever the species was observed, even though 
    # some of those observations may be false positives. 
    z.inits <- apply(y.big, 1, max, na.rm = TRUE)
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # beta -----------------------
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p.occ & length(beta.inits) != 1) {
      if (p.occ == 1) {
        stop(paste("initial values for beta must be of length ", p.occ,
                   sep = ""))
        
      } else {
        stop(paste("initial values for beta must be of length ", p.occ, " or 1",
                   sep = ""))
      }
    }
    if (length(beta.inits) != p.occ) {
      beta.inits <- rep(beta.inits, p.occ)
    }
  } else {
    beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
    if (verbose) {
      message('beta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n')
    }
  }
  # gamma -----------------------
  if ("gamma" %in% names(inits)) {
    gamma.inits <- inits[["gamma"]]
    if (length(gamma.inits) != p.tp.det & length(gamma.inits) != 1) {
      if (p.tp.det == 1) {
        stop(paste("initial values for gamma must be of length ", p.tp.det,
                   sep = ""))
      } else {
        stop(paste("initial values for gamma must be of length ", p.tp.det, " or 1",
                   sep = ""))
      }
    }
    if (length(gamma.inits) != p.tp.det) {
      gamma.inits <- rep(gamma.inits, p.tp.det)
    }
  } else {
    gamma.inits <- rnorm(p.tp.det, mu.gamma, sqrt(sigma.gamma))
    if (verbose) {
      message("gamma is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  # delta -----------------------
  if ("delta" %in% names(inits)) {
    delta.inits <- inits[["delta"]]
    if (length(delta.inits) != p.fp.det & length(delta.inits) != 1) {
      if (p.fp.det == 1) {
        stop(paste("initial values for delta must be of length ", p.fp.det,
                   sep = ""))
      } else {
        stop(paste("initial values for delta must be of length ", p.fp.det, " or 1",
                   sep = ""))
      }
    }
    if (length(delta.inits) != p.fp.det) {
      delta.inits <- rep(delta.inits, p.fp.det)
    }
  } else {
    delta.inits <- rnorm(p.fp.det, mu.delta, sqrt(sigma.delta))
    if (verbose) {
      message("delta is not specified in initial values.\nSetting initial values to random values from the prior distribution\n")
    }
  }
  
  # sigma.sq.psi -------------------
  if (p.occ.re > 0) {
    if ("sigma.sq.psi" %in% names(inits)) {
      sigma.sq.psi.inits <- inits[["sigma.sq.psi"]]
      if (length(sigma.sq.psi.inits) != p.occ.re & length(sigma.sq.psi.inits) != 1) {
        if (p.occ.re == 1) {
          stop(paste("initial values for sigma.sq.psi must be of length ", p.occ.re, 
                     sep = ""))
        } else {
          stop(paste("initial values for sigma.sq.psi must be of length ", p.occ.re, 
                     " or 1", sep = ""))
        }
      }
      if (length(sigma.sq.psi.inits) != p.occ.re) {
        sigma.sq.psi.inits <- rep(sigma.sq.psi.inits, p.occ.re)
      }
    } else {
      sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
      if (verbose) {
        message("sigma.sq.psi is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    beta.star.indx <- rep(0:(p.occ.re - 1), n.occ.re.long)
    beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
  } else {
    sigma.sq.psi.inits <- 0
    beta.star.indx <- 0
    beta.star.inits <- 0
  }
  # sigma.sq.mu ------------------
  if (p.tp.det.re > 0) {
    if ("sigma.sq.mu" %in% names(inits)) {
      sigma.sq.mu.inits <- inits[["sigma.sq.mu"]]
      if (length(sigma.sq.mu.inits) != p.tp.det.re & length(sigma.sq.mu.inits) != 1) {
        if (p.tp.det.re == 1) {
          stop(paste("initial values for sigma.sq.mu must be of length ", p.tp.det.re, 
                     sep = ""))
        } else {
          stop(paste("initial values for sigma.sq.mu must be of length ", p.tp.det.re, 
                     " or 1", sep = ""))
          
        }
      }
      if (length(sigma.sq.mu.inits) != p.tp.det.re) {
        sigma.sq.mu.inits <- rep(sigma.sq.mu.inits, p.tp.det.re)
      }
    } else {
      sigma.sq.mu.inits <- runif(p.tp.det.re, 0.5, 4)
      if (verbose) {
        message("sigma.sq.mu is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    gamma.star.indx <- rep(0:(p.tp.det.re - 1), n.tp.det.re.long)
    gamma.star.inits <- rnorm(n.tp.det.re, 0, sqrt(sigma.sq.mu.inits[gamma.star.indx + 1]))
  } else {
    sigma.sq.mu.inits <- 0
    gamma.star.indx <- 0
    gamma.star.inits <- 0
  }
  # sigma.sq.omega ------------------
  if (p.fp.det.re > 0) {
    if ("sigma.sq.omega" %in% names(inits)) {
      sigma.sq.omega.inits <- inits[["sigma.sq.omega"]]
      if (length(sigma.sq.omega.inits) != p.fp.det.re & length(sigma.sq.omega.inits) != 1) {
        if (p.fp.det.re == 1) {
          stop(paste("initial values for sigma.sq.omega must be of length ", p.fp.det.re, 
                     sep = ""))
        } else {
          stop(paste("initial values for sigma.sq.omega must be of length ", p.fp.det.re, 
                     " or 1", sep = ""))
          
        }
      }
      if (length(sigma.sq.omega.inits) != p.fp.det.re) {
        sigma.sq.omega.inits <- rep(sigma.sq.omega.inits, p.fp.det.re)
      }
    } else {
      sigma.sq.omega.inits <- runif(p.fp.det.re, 0.5, 4)
      if (verbose) {
        message("sigma.sq.omega is not specified in initial values.\nSetting initial values to random values between 0.5 and 10\n")
      }
    }
    delta.star.indx <- rep(0:(p.fp.det.re - 1), n.fp.det.re.long)
    delta.star.inits <- rnorm(n.fp.det.re, 0, sqrt(sigma.sq.omega.inits[delta.star.indx + 1]))
  } else {
    sigma.sq.omega.inits <- 0
    delta.star.indx <- 0
    delta.star.inits <- 0
  }
  # D ---------------------------------
  # NOTE: you look for d since you call tolower() on the user-supplied inits list. 
  if ("d" %in% names(inits)) {
    D.inits <- inits$d
    if (!is.matrix(D.inits)) {
      stop(paste("initial values for D must be a matrix with", J, "rows and", 
                 K.max, "columns"))
    }
    if (nrow(D.inits) != J | ncol(D.inits) != K.max) {
      stop(paste("initial values for D must be a matrix with", J, "rows and", 
                 K.max, "columns"))
    }
    D.inits <- D.inits[!is.na(d.big)]
  } else {
    # Default is set to the number of true positives in the validation data set
    D.inits <- d.big
    if (verbose) {
      message("D is not specified in initial values.\nSetting initial values based on observed validation data\n")
    }
    D.inits <- D.inits[!is.na(d.big)]
  }
  
  # Should initial values be fixed --
  if ("fix" %in% names(inits)) {
    fix.inits <- inits[["fix"]]
    if ((fix.inits != TRUE) & (fix.inits != FALSE)) {
      stop(paste("inits$fix must take value TRUE or FALSE"))
    }
  } else {
    fix.inits <- FALSE
  }
  if (verbose & fix.inits & (n.chains > 1)) {
    message("Fixing initial values across all chains\n")
  }
 
  # Tuning values --------------------------------------------------------- 
  if (missing(tuning)) {
    gamma.tuning <- rep(1, p.tp.det)
    delta.tuning <- rep(1, p.fp.det)
    gamma.star.tuning <- rep(1, n.tp.det.re)
    delta.star.tuning <- rep(1, n.fp.det.re)
  } else {
    names(tuning) <- tolower(names(tuning))
    # gamma ---------------------------
    if (!"gamma" %in% names(tuning)) {
      stop("gamma must be specified in tuning value list")
    }
    gamma.tuning <- tuning$gamma
    if (length(gamma.tuning) != 1 & length(gamma.tuning) != p.tp.det) {
      stop(paste("gamma tuning must be a single value or a vector of length ",
                 p.tp.det, sep = ''))
    }
    if (length(gamma.tuning) == 1) {
      gamma.tuning <- rep(gamma.tuning, p.tp.det)
    }
    # delta ---------------------------
    if (!"delta" %in% names(tuning)) {
      stop("delta must be specified in tuning value list")
    }
    delta.tuning <- tuning$delta
    if (length(delta.tuning) != 1 & length(delta.tuning) != p.fp.det) {
      stop(paste("delta tuning must be a single value or a vector of length ",
                 p.fp.det, sep = ''))
    }
    if (length(delta.tuning) == 1) {
      delta.tuning <- rep(delta.tuning, p.fp.det)
    }
    if (p.tp.det.re > 0) {
      # gamma.star ---------------------------
      if(!"gamma.star" %in% names(tuning)) {
        stop("gamma.star must be specified in tuning value list")
      }
      gamma.star.tuning <- tuning$gamma.star
      if (length(gamma.star.tuning) != 1) {
        stop("gamma.star tuning must be a single value")
      }
      gamma.star.tuning <- rep(gamma.star.tuning, n.tp.det.re)
    } else {
      gamma.star.tuning <- NULL
    }
    if (p.fp.det.re > 0) {
      # delta.star ---------------------------
      if(!"delta.star" %in% names(tuning)) {
        stop("delta.star must be specified in tuning value list")
      }
      delta.star.tuning <- tuning$delta.star
      if (length(delta.star.tuning) != 1) {
        stop("delta.star tuning must be a single value")
      }
      delta.star.tuning <- rep(delta.star.tuning, n.fp.det.re)
    } else {
      delta.star.tuning <- NULL
    }
  }
  tuning.c <- log(c(gamma.tuning, delta.tuning, gamma.star.tuning,
                    delta.star.tuning))
  
  
  curr.chain <- 1
  # Set storage for all variables ---------------------------------------
  storage.mode(y) <- "double"
  storage.mode(z.inits) <- "double"
  storage.mode(n.val) <- "double"
  storage.mode(d) <- "double"
  storage.mode(X.mu) <- "double"
  storage.mode(X.omega) <- "double"
  storage.mode(X) <- "double"
  consts <- c(J, n.obs, p.occ, p.occ.re, n.occ.re, p.tp.det, p.tp.det.re, 
              n.tp.det.re, p.fp.det, p.fp.det.re, n.fp.det.re)
  storage.mode(consts) <- "integer"
  storage.mode(K) <- "double"
  storage.mode(beta.inits) <- "double"
  storage.mode(gamma.inits) <- "double"
  storage.mode(delta.inits) <- "double"
  storage.mode(D.inits) <- "double"
  storage.mode(z.long.indx) <- "integer"
  storage.mode(mu.beta) <- "double"
  storage.mode(Sigma.beta) <- "double"
  storage.mode(mu.gamma) <- "double"
  storage.mode(Sigma.gamma) <- "double"
  storage.mode(mu.delta) <- "double"
  storage.mode(Sigma.delta) <- "double"
  storage.mode(tuning.c) <- "double"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"
  storage.mode(n.omp.threads) <- "integer"
  storage.mode(verbose) <- "integer"
  storage.mode(n.report) <- "integer"
  chain.info <- c(curr.chain, n.chains)
  storage.mode(chain.info) <- "integer"
  n.post.samples <- length(seq(from = n.burn + 1, 
                               to = n.samples, 
                               by = as.integer(n.thin)))
  storage.mode(n.post.samples) <- "integer"
  samples.info <- c(n.burn, n.thin, n.post.samples)
  storage.mode(samples.info) <- "integer"
  # For true positive detection random effects
  storage.mode(X.mu.re) <- "integer"
  gamma.level.indx <- sort(unique(c(X.mu.re)))
  storage.mode(gamma.level.indx) <- "integer"
  storage.mode(n.tp.det.re.long) <- "integer"
  storage.mode(sigma.sq.mu.inits) <- "double"
  storage.mode(sigma.sq.mu.a) <- "double"
  storage.mode(sigma.sq.mu.b) <- "double"
  storage.mode(gamma.star.inits) <- "double"
  storage.mode(gamma.star.indx) <- "integer"
  # For false positive detection random effects
  storage.mode(X.omega.re) <- "integer"
  delta.level.indx <- sort(unique(c(X.omega.re)))
  storage.mode(delta.level.indx) <- "integer"
  storage.mode(n.fp.det.re.long) <- "integer"
  storage.mode(sigma.sq.omega.inits) <- "double"
  storage.mode(sigma.sq.omega.a) <- "double"
  storage.mode(sigma.sq.omega.b) <- "double"
  storage.mode(delta.star.inits) <- "double"
  storage.mode(delta.star.indx) <- "integer"
  # For occurrence random effects
  storage.mode(X.re) <- "integer"
  beta.level.indx <- sort(unique(c(X.re)))
  storage.mode(beta.level.indx) <- "integer"
  storage.mode(sigma.sq.psi.inits) <- "double"
  storage.mode(sigma.sq.psi.a) <- "double"
  storage.mode(sigma.sq.psi.b) <- "double"
  storage.mode(beta.star.inits) <- "double"
  storage.mode(beta.star.indx) <- "integer"
  
  # Fit the model -------------------------------------------------------
  out.tmp <- list()
  out <- list()
  for (i in 1:n.chains) {
    # Change initial values if i > 1
    if ((i > 1) & (!fix.inits)) {
      beta.inits <- rnorm(p.occ, mu.beta, sqrt(sigma.beta))
      gamma.inits <- rnorm(p.tp.det, mu.gamma, sqrt(sigma.gamma))
      delta.inits <- rnorm(p.fp.det, mu.delta, sqrt(sigma.delta))
      if (p.occ.re > 0) {
        sigma.sq.psi.inits <- runif(p.occ.re, 0.5, 10)
        beta.star.inits <- rnorm(n.occ.re, 0, sqrt(sigma.sq.psi.inits[beta.star.indx + 1]))
      }
      if (p.tp.det.re > 0) {
        sigma.sq.mu.inits <- runif(p.tp.det.re, 0.5, 3)
        gamma.star.inits <- rnorm(n.tp.det.re, 0, sqrt(sigma.sq.mu.inits[gamma.star.indx + 1]))
      }
      if (p.fp.det.re > 0) {
        sigma.sq.omega.inits <- runif(p.fp.det.re, 0.5, 3)
        delta.star.inits <- rnorm(n.fp.det.re, 0, sqrt(sigma.sq.omega.inits[delta.star.indx + 1]))
      }
    }
    storage.mode(chain.info) <- "integer"
    # Run the model in C
    out.tmp[[i]] <- .Call("fpPGOcc", y, X, X.mu, X.omega, X.re, X.mu.re, X.omega.re, 
                          d, n.val, consts, K, n.occ.re.long, n.tp.det.re.long, 
                          n.fp.det.re.long, beta.inits, gamma.inits, delta.inits,
                          sigma.sq.psi.inits, sigma.sq.mu.inits, sigma.sq.omega.inits, 
                          beta.star.inits, gamma.star.inits, delta.star.inits, z.inits, 
                          D.inits, z.long.indx, beta.star.indx, beta.level.indx, 
                          gamma.star.indx, gamma.level.indx, delta.star.indx, 
                          delta.level.indx, mu.beta, mu.gamma, mu.delta, 
                          Sigma.beta, Sigma.gamma, Sigma.delta, sigma.sq.psi.a, 
                          sigma.sq.psi.b, sigma.sq.mu.a, sigma.sq.mu.b, 
                          sigma.sq.omega.a, sigma.sq.omega.b, tuning.c, n.batch,
                          batch.length, accept.rate, n.omp.threads, verbose, 
                          n.report, samples.info, chain.info)
    chain.info[1] <- chain.info[1] + 1
  } # i   
  # Calculate R-Hat ---------------
  out <- list()
  out$rhat <- list()
  if (n.chains > 1) {
    # as.vector removes the "Upper CI" when there is only 1 variable. 
    out$rhat$beta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      mcmc(t(a$beta.samples)))), 
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$gamma <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      mcmc(t(a$gamma.samples)))), 
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    out$rhat$delta <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
      mcmc(t(a$delta.samples)))), 
      autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    if (p.tp.det.re > 0) {
      out$rhat$sigma.sq.mu <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        mcmc(t(a$sigma.sq.mu.samples)))), 
        autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
    if (p.fp.det.re > 0) {
      out$rhat$sigma.sq.omega <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        mcmc(t(a$sigma.sq.omega.samples)))), 
        autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
    if (p.occ.re > 0) {
      out$rhat$sigma.sq.psi <- as.vector(gelman.diag(mcmc.list(lapply(out.tmp, function(a) 
        mcmc(t(a$sigma.sq.psi.samples)))), 
        autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])
    }
  } else {
    out$rhat$beta <- rep(NA, p.occ)
    out$rhat$gamma <- rep(NA, p.tp.det)
    out$rhat$delta <- rep(NA, p.fp.det)
    if (p.tp.det.re > 0) {
      out$rhat$sigma.sq.mu <- rep(NA, p.tp.det.re)
    }
    if (p.fp.det.re > 0) {
      out$rhat$sigma.sq.omega <- rep(NA, p.fp.det.re)
    }
    if (p.occ.re > 0) {
      out$rhat$sigma.sq.psi <- rep(NA, p.occ.re)
    }
  }
  # Put everything into MCMC objects
  out$beta.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$beta.samples))))
  colnames(out$beta.samples) <- x.names
  out$gamma.samples <- mcmc(do.call(rbind, 
                                    lapply(out.tmp, function(a) t(a$gamma.samples))))
  colnames(out$gamma.samples) <- x.mu.names
  out$delta.samples <- mcmc(do.call(rbind, 
                                    lapply(out.tmp, function(a) t(a$delta.samples))))
  colnames(out$delta.samples) <- x.omega.names
  out$z.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$z.samples))))
  # TODO: this needs to be fixed eventually to give the observation values in the right 
  #       array shape. 
  out$D.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$D.samples))))
  out$psi.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$psi.samples))))
  out$like.samples <- mcmc(do.call(rbind, lapply(out.tmp, function(a) t(a$like.samples))))
  if (p.occ.re > 0) {
    out$sigma.sq.psi.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.psi.samples))))
    colnames(out$sigma.sq.psi.samples) <- x.re.names
    out$beta.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$beta.star.samples))))
    tmp.names <- unlist(re.level.names)
    beta.star.names <- paste(rep(x.re.names, n.occ.re.long), tmp.names, sep = '-')
    colnames(out$beta.star.samples) <- beta.star.names
    out$re.level.names <- re.level.names
  }
  if (p.tp.det.re > 0) {
    out$sigma.sq.mu.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.mu.samples))))
    colnames(out$sigma.sq.mu.samples) <- x.mu.re.names
    out$gamma.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$gamma.star.samples))))
    tmp.names <- unlist(mu.re.level.names)
    gamma.star.names <- paste(rep(x.mu.re.names, n.tp.det.re.long), tmp.names, sep = '-')
    colnames(out$gamma.star.samples) <- gamma.star.names
    out$mu.re.level.names <- mu.re.level.names
  }
  if (p.fp.det.re > 0) {
    out$sigma.sq.omega.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$sigma.sq.omega.samples))))
    colnames(out$sigma.sq.omega.samples) <- x.omega.re.names
    out$delta.star.samples <- mcmc(
      do.call(rbind, lapply(out.tmp, function(a) t(a$delta.star.samples))))
    tmp.names <- unlist(omega.re.level.names)
    delta.star.names <- paste(rep(x.omega.re.names, n.fp.det.re.long), tmp.names, sep = '-')
    colnames(out$delta.star.samples) <- delta.star.names
    out$omega.re.level.names <- omega.re.level.names
  }
  # Calculate effective sample sizes
  out$ESS <- list()
  out$ESS$beta <- effectiveSize(out$beta.samples)
  out$ESS$gamma <- effectiveSize(out$gamma.samples)
  if (p.tp.det.re > 0) {
    out$ESS$sigma.sq.mu <- effectiveSize(out$sigma.sq.mu.samples)
  }
  if (p.fp.det.re > 0) {
    out$ESS$sigma.sq.omega <- effectiveSize(out$sigma.sq.omega.samples)
  }
  if (p.occ.re > 0) {
    out$ESS$sigma.sq.psi <- effectiveSize(out$sigma.sq.psi.samples)
  }
  out$X <- X
  out$X.mu <- X.mu
  out$X.omega <- X.omega
  out$X.re <- X.re
  out$X.mu.re <- X.mu.re
  out$X.omega.re <- X.omega.re
  out$y <- y.big
  out$n.samples <- batch.length * n.batch
  out$call <- cl
  out$n.post <- n.post.samples
  out$n.thin <- n.thin
  out$n.burn <- n.burn
  out$n.chains <- n.chains
  if (p.tp.det.re > 0) {
    out$muRE <- TRUE
  } else {
    out$muRE <- FALSE
  }
  if (p.fp.det.re > 0) {
    out$omegaRE <- TRUE
  } else {
    out$omegaRE <- FALSE
  }
  if (p.occ.re > 0) {
    out$psiRE <- TRUE
  } else {
    out$psiRE <- FALSE
  }
  
  return(out)
}
