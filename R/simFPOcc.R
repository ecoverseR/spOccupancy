simFPOcc <- function(J.x, J.y, n.rep, n.rep.max, beta, gamma, delta, 
                     psi.RE = list(), mu.RE = list(), omega.RE = omega.RE, 
                     sp = FALSE, svc.cols = 1, cov.model, sigma.sq, phi, nu, 
                     x.positive = FALSE, grid, prop.val, kappa, family = "Poisson", ...) {

  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }

  # Subroutines -----------------------------------------------------------
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }

  # Check function inputs -------------------------------------------------
  # J.x -------------------------------
  if (missing(J.x)) {
    stop("J.x must be specified")
  }
  if (length(J.x) != 1) {
    stop("J.x must be a single numeric value.")
  }
  # J.y -------------------------------
  if (missing(J.y)) {
    stop("J.y must be specified")
  }
  if (length(J.y) != 1) {
    stop("J.y must be a single numeric value.")
  }
  J <- J.x * J.y
  # n.rep -----------------------------
  if (missing(n.rep)) {
    stop("n.rep must be specified.")
  }
  if (length(n.rep) != J) {
    stop(paste("n.rep must be a vector of length ", J, sep = ''))
  }
  if (missing(n.rep.max)) {
    n.rep.max <- max(n.rep)
  }
  # beta ------------------------------
  if (missing(beta)) {
    stop("beta must be specified.")
  }
  # gamma -----------------------------
  if (missing(gamma)) {
    stop("gamma must be specified.")
  }
  # delta -----------------------------
  if (missing(delta)) {
    stop("delta must be specified.")
  }
  # kappa -----------------------------
  if (family == 'NB') {
    if (missing(kappa)) {
      stop("kappa (overdispersion parameter) must be specified when family = 'NB'.")
    }
  }
  if (family == 'Poisson' & !missing(kappa)) {
    message("overdispersion parameter (kappa) is ignored when family == 'Poisson'")
  }
  # psi.RE ----------------------------
  names(psi.RE) <- tolower(names(psi.RE))
  if (!is.list(psi.RE)) {
    stop("if specified, psi.RE must be a list with tags 'levels' and 'sigma.sq.psi'")
  }
  if (length(names(psi.RE)) > 0) {
    if (!'sigma.sq.psi' %in% names(psi.RE)) {
      stop("sigma.sq.psi must be a tag in psi.RE with values for the occurrence random effect variances")
    }
    if (!'levels' %in% names(psi.RE)) {
      stop("levels must be a tag in psi.RE with the number of random effect levels for each occurrence random intercept.")
    }
    if (!'beta.indx' %in% names(psi.RE)) {
      psi.RE$beta.indx <- list()
      for (i in 1:length(psi.RE$sigma.sq.psi)) {
        psi.RE$beta.indx[[i]] <- 1
      }
    }
  }
  # mu.RE ----------------------------
  names(mu.RE) <- tolower(names(mu.RE))
  if (!is.list(mu.RE)) {
    stop("if specified, mu.RE must be a list with tags 'levels' and 'sigma.sq.mu'")
  }
  if (length(names(mu.RE)) > 0) {
    if (!'sigma.sq.mu' %in% names(mu.RE)) {
      stop("sigma.sq.mu must be a tag in mu.RE with values for the true positive detection random effect variances")
    }
    if (!'levels' %in% names(mu.RE)) {
      stop("levels must be a tag in mu.RE with the number of random effect levels for each true positive detection rate random intercept.")
    }
    if (!'gamma.indx' %in% names(mu.RE)) {
      mu.RE$gamma.indx <- list()
      for (i in 1:length(mu.RE$sigma.sq.mu)) {
        mu.RE$gamma.indx[[i]] <- 1
      }
    }
  }
  # omega.RE --------------------------
  names(omega.RE) <- tolower(names(omega.RE))
  if (!is.list(omega.RE)) {
    stop("if specified, omega.RE must be a list with tags 'levels' and 'sigma.sq.omega'")
  }
  if (length(names(omega.RE)) > 0) {
    if (!'sigma.sq.omega' %in% names(omega.RE)) {
      stop("sigma.sq.omega must be a tag in omega.RE with values for the false positive detection random effect variances")
    }
    if (!'levels' %in% names(omega.RE)) {
      stop("levels must be a tag in omega.RE with the number of random effect levels for each false positive detection rate random intercept.")
    }
    if (!'delta.indx' %in% names(omega.RE)) {
      omega.RE$delta.indx <- list()
      for (i in 1:length(omega.RE$sigma.sq.omega)) {
        omega.RE$delta.indx[[i]] <- 1
      }
    }
  }
  if (length(svc.cols) > 1 & !sp) {
    stop("if simulating data with spatially-varying coefficients, set sp = TRUE")
  }
  # Spatial parameters ----------------
  p.svc <- length(svc.cols)
  if (sp) {
    if(missing(sigma.sq)) {
      stop("sigma.sq must be specified when sp = TRUE")
    }
    if(missing(phi)) {
      stop("phi must be specified when sp = TRUE")
    }
    if(missing(cov.model)) {
      stop("cov.model must be specified when sp = TRUE")
    }
    cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
    if(! cov.model %in% cov.model.names){
      stop("specified cov.model '",cov.model,"' is not a valid option; choose from ", 
           paste(cov.model.names, collapse=", ", sep="") ,".")
    }
    if (cov.model == 'matern' & missing(nu)) {
      stop("nu must be specified when cov.model = 'matern'")
    }
    if (length(phi) != p.svc) {
      stop("phi must have the same number of elements as svc.cols")
    }
    if (length(sigma.sq) != p.svc) {
      stop("sigma.sq must have the same number of elements as svc.cols")
    }
    if (cov.model == 'matern') {
      if (length(nu) != p.svc) {
        stop("nu must have the same number of elements as svc.cols")
      }
    }
  }
  # Grid for spatial REs that doesn't match the sites ---------------------
  if (!missing(grid) & sp) {
    if (!is.atomic(grid)) {
      stop("grid must be a vector")
    }
    if (length(grid) != J) {
      stop(paste0("grid must be of length ", J))
    }
  } else {
    grid <- 1:J
  }

  # Subroutines -----------------------------------------------------------
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

  # Matrix of spatial locations
  s.x <- seq(0, 1, length.out = J.x)
  s.y <- seq(0, 1, length.out = J.y)
  coords.full <- as.matrix(expand.grid(s.x, s.y))
  coords <- cbind(tapply(coords.full[, 1], grid, mean), 
                  tapply(coords.full[, 2], grid, mean)) 

  # Form occupancy covariates (if any) ------------------------------------
  n.beta <- length(beta)
  X <- matrix(1, nrow = J, ncol = n.beta)
  if (n.beta > 1) {
    for (i in 2:n.beta) {
      X[, i] <- rnorm(J)
    } # i
    if (x.positive) {
      for (i in 2:n.beta) {
        X[, i] <- runif(J, 0, 1)
      }
    }
  }

  # Form detection covariates (if any) ------------------------------------
  # Get index of surveyed replicates for each site. 
  rep.indx <- list()
  for (j in 1:J) {
    rep.indx[[j]] <- sample(1:n.rep.max, n.rep[j], replace = FALSE)
  }
  # True positive rate ----------------
  n.gamma <- length(gamma)
  X.mu <- array(NA, dim = c(J, n.rep.max, n.gamma))
  X.mu[, , 1] <- 1
  if (n.gamma > 1) {
    for (i in 2:n.gamma) {
      for (j in 1:J) {
        X.mu[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
      } # j
    } # i
  }
  
  # False positive rate ----------------
  n.delta <- length(delta)
  X.omega <- array(NA, dim = c(J, n.rep.max, n.delta))
  X.omega[, , 1] <- 1
  if (n.delta > 1) {
    for (i in 2:n.delta) {
      for (j in 1:J) {
        X.omega[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
      } # j
    } # i
  }
  
  # Simulate spatial random effect ----------------------------------------
  # Matrix of spatial locations
  if (sp) {
    J <- nrow(coords.full)
    w.mat.full <- matrix(NA, J, p.svc)
    if (cov.model == 'matern') {
      theta <- cbind(phi, nu)
    } else {
      theta <- as.matrix(phi)
    }
    w.mat <- matrix(NA, nrow(coords), p.svc)
    for (i in 1:p.svc) {
      Sigma.full <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
      w.mat[, i] <- rmvn(1, rep(0, nrow(Sigma.full)), Sigma.full)
      # Random spatial process
      w.mat.full[, i] <- w.mat[grid, i]
    }
    X.w <- X[, svc.cols, drop = FALSE]
    # Convert w to a J*ptilde x 1 vector, sorted so that the p.svc values for 
    # each site are given, then the next site, then the next, etc.
    w <- c(t(w.mat.full))
    # Create X.tilde, which is a J x J*p.tilde matrix. 
    X.tilde <- matrix(0, J, J * p.svc)
    # Fill in the matrix
    for (j in 1:J) {
      X.tilde[j, ((j - 1) * p.svc + 1):(j * p.svc)] <- X.w[j, ]
    }
  } else {
    w.mat.full <- NA
    w.mat <- NA
    X.w <- NA
  }

  # Random effects --------------------------------------------------------
  if (length(psi.RE) > 0) {
    p.occ.re <- length(unlist(psi.RE$beta.indx))
    tmp <- sapply(psi.RE$beta.indx, length)
    re.col.indx <- unlist(lapply(1:length(psi.RE$beta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.psi <- psi.RE$sigma.sq.psi[re.col.indx]
    n.occ.re.long <- psi.RE$levels[re.col.indx]
    n.occ.re <- sum(n.occ.re.long)
    beta.star.indx <- rep(1:p.occ.re, n.occ.re.long)
    beta.star <- rep(0, n.occ.re)
    X.random <- X[, unlist(psi.RE$beta.indx), drop = FALSE]
    n.random <- ncol(X.random)
    X.re <- matrix(NA, J, length(psi.RE$levels))
    for (i in 1:length(psi.RE$levels)) {
      X.re[, i] <- sample(1:psi.RE$levels[i], J, replace = TRUE)  
    }
    indx.mat <- X.re[, re.col.indx, drop = FALSE]
    for (i in 1:p.occ.re) {
      beta.star[which(beta.star.indx == i)] <- rnorm(n.occ.re.long[i], 0, 
						     sqrt(sigma.sq.psi[i]))
    }
    if (length(psi.RE$levels) > 1) {
      for (j in 2:length(psi.RE$levels)) {
        X.re[, j] <- X.re[, j] + max(X.re[, j - 1], na.rm = TRUE)
      }
    }
    if (p.occ.re > 1) {
      for (j in 2:p.occ.re) {
        indx.mat[, j] <- indx.mat[, j] + max(indx.mat[, j - 1], na.rm = TRUE)
      }
    }
    beta.star.sites <- rep(NA, J)
    for (j in 1:J) {
      beta.star.sites[j] <- beta.star[indx.mat[j, , drop = FALSE]] %*% t(X.random[j, , drop = FALSE])
    }
  } else {
    X.re <- NA
    beta.star <- NA
  }
  # False positive rate random effects
  if (length(omega.RE) > 0) {
    p.fp.re <- length(unlist(omega.RE$delta.indx))
    tmp <- sapply(omega.RE$delta.indx, length)
    omega.re.col.indx <- unlist(lapply(1:length(omega.RE$delta.indx), function(a) rep(a, tmp[a])))
    sigma.sq.omega <- omega.RE$sigma.sq.omega[omega.re.col.indx]
    n.fp.re.long <- omega.RE$levels[omega.re.col.indx]
    n.fp.re <- sum(n.fp.re.long)
    delta.star.indx <- rep(1:p.fp.re, n.fp.re.long)
    delta.star <- rep(0, n.fp.re)
    X.omega.random <- X.omega[, , unlist(omega.RE$delta.indx), drop = FALSE]
    X.omega.re <- array(NA, dim = c(J, n.rep.max, length(omega.RE$levels)))
    for (i in 1:length(omega.RE$levels)) {
      X.omega.re[, , i] <- matrix(sample(1:omega.RE$levels[i], J * n.rep.max, replace = TRUE), 
		              J, n.rep.max)	      
    }
    for (i in 1:p.fp.re) {
      delta.star[which(delta.star.indx == i)] <- rnorm(n.fp.re.long[i], 0, sqrt(sigma.sq.omega[i]))
    }
    X.omega.re <- X.omega.re[, , omega.re.col.indx, drop = FALSE]
    for (j in 1:J) {
      X.omega.re[j, -rep.indx[[j]], ] <- NA
    }
    if (p.fp.re > 1) {
      for (j in 2:p.fp.re) {
        X.omega.re[, , j] <- X.omega.re[, , j] + max(X.omega.re[, , j - 1], na.rm = TRUE) 
      }
    }
    delta.star.sites <- matrix(NA, J, n.rep.max)
    for (j in 1:J) {
      for (k in rep.indx[[j]]) {
        delta.star.sites[j, k] <- delta.star[X.omega.re[j, k, ]] %*% X.omega.random[j, k, ]
      }
    }
  } else {
    X.omega.re <- NA
    delta.star <- NA
  }
  # True positive rate random effects
  if (length(mu.RE) > 0) {
    p.tp.re <- length(unlist(mu.RE$gamma.indx))
    tmp <- sapply(mu.RE$gamma.indx, length)
    mu.re.col.indx <- unlist(lapply(1:length(mu.RE$gamma.indx), function(a) rep(a, tmp[a])))
    sigma.sq.mu <- mu.RE$sigma.sq.mu[mu.re.col.indx]
    n.tp.re.long <- mu.RE$levels[mu.re.col.indx]
    n.tp.re <- sum(n.tp.re.long)
    gamma.star.indx <- rep(1:p.tp.re, n.tp.re.long)
    gamma.star <- rep(0, n.tp.re)
    X.mu.random <- X.mu[, , unlist(mu.RE$gamma.indx), drop = FALSE]
    X.mu.re <- array(NA, dim = c(J, n.rep.max, length(mu.RE$levels)))
    for (i in 1:length(mu.RE$levels)) {
      X.mu.re[, , i] <- matrix(sample(1:mu.RE$levels[i], J * n.rep.max, replace = TRUE), 
		              J, n.rep.max)	      
    }
    for (i in 1:p.tp.re) {
      gamma.star[which(gamma.star.indx == i)] <- rnorm(n.tp.re.long[i], 0, sqrt(sigma.sq.mu[i]))
    }
    X.mu.re <- X.mu.re[, , mu.re.col.indx, drop = FALSE]
    for (j in 1:J) {
      X.mu.re[j, -rep.indx[[j]], ] <- NA
    }
    if (p.tp.re > 1) {
      for (j in 2:p.tp.re) {
        X.mu.re[, , j] <- X.mu.re[, , j] + max(X.mu.re[, , j - 1], na.rm = TRUE) 
      }
    }
    gamma.star.sites <- matrix(NA, J, n.rep.max)
    for (j in 1:J) {
      for (k in rep.indx[[j]]) {
        gamma.star.sites[j, k] <- gamma.star[X.mu.re[j, k, ]] %*% X.mu.random[j, k, ]
      }
    }
  } else {
    X.mu.re <- NA
    gamma.star <- NA
  }
  # Latent Occupancy Process ----------------------------------------------
  if (sp) {
    if (length(psi.RE) > 0) {
      psi <- logit.inv(X %*% as.matrix(beta) + X.tilde %*% w + beta.star.sites)
    } else {
      psi <- logit.inv(X %*% as.matrix(beta) + X.tilde %*% w)
    }
  } else {
    if (length(psi.RE) > 0) {
      psi <- logit.inv(X %*% as.matrix(beta) + beta.star.sites)
    } else {
      psi <- logit.inv(X %*% as.matrix(beta))
    }
  }
  z <- rbinom(J, 1, psi)

  # Data Formation --------------------------------------------------------
  # Observed count data 
  y <- matrix(NA, nrow = J, ncol = n.rep.max)
  # Latent true positive rate. 
  mu <- matrix(NA, nrow = J, ncol = n.rep.max)
  # Latent false positive rate.
  omega <- matrix(NA, nrow = J, ncol = n.rep.max)
  # Number of validated samples (observed)
  n.val <- matrix(NA, nrow = J, ncol = n.rep.max)
  # Number of true positive of validated samples (observed)
  d <- matrix(NA, nrow = J, ncol = n.rep.max)
  # Latent number of true positives
  D <- matrix(NA, nrow = J, ncol = n.rep.max)
  # Latent number of false positives
  Q <- matrix(NA, nrow = J, ncol = n.rep.max)

  # Generate data and latent vocalization processes
  for (j in 1:J) {
    # TP rate
    if (length(mu.RE) > 0) {
      mu[j, rep.indx[[j]]] <- exp(X.mu[j, rep.indx[[j]], ] %*% as.matrix(gamma) + 
				    gamma.star.sites[j, rep.indx[[j]]])
    } else {
      mu[j, rep.indx[[j]]] <- exp(X.mu[j, rep.indx[[j]], ] %*% as.matrix(gamma))
    }
    # FP rate
    if (length(omega.RE) > 0) {
      omega[j, rep.indx[[j]]] <- exp(X.omega[j, rep.indx[[j]], ] %*% as.matrix(delta) + 
				    delta.star.sites[j, rep.indx[[j]]])
    } else {
      omega[j, rep.indx[[j]]] <- exp(X.omega[j, rep.indx[[j]], ] %*% as.matrix(delta))
    }
    # Get total vocalizations/images at a site
    if (family == "Poisson") {
      y[j, rep.indx[[j]]] <- rpois(n.rep[j], mu[j, rep.indx[[j]]] * z[j] + omega[j, rep.indx[[j]]])
    } else if (family == "NB") {
      y[j, rep.indx[[j]]] <- rnbinom(n.rep[j], size = kappa, mu = mu[j, rep.indx[[j]]] * z[j] + omega[j, rep.indx[[j]]])
    }
    # Get total true positives at each site. 
    D[j, rep.indx[[j]]] <- rbinom(n.rep[j], 
                                  y[j, rep.indx[[j]]],
                                  (mu[j, rep.indx[[j]]] * z[j]) / (mu[j, rep.indx[[j]]] * z[j] + omega[j, rep.indx[[j]]])) 

    # Get total false positives
    Q[j, rep.indx[[j]]] <- y[j, rep.indx[[j]]] - D[j, rep.indx[[j]]]
    
    # Number of vocalizations at the current site/rep to be validated
    n.val[j, rep.indx[[j]]] <- round(prop.val * y[j, rep.indx[[j]]])
    d[j, rep.indx[[j]]] <- rhyper(n.rep[j], D[j, rep.indx[[j]]], Q[j, rep.indx[[j]]], n.val[j, rep.indx[[j]]])
  }


  return(
    list(X = X, X.mu = X.mu, X.omega = X.omega, coords = coords, coords.full = coords.full,
         w = w.mat.full, w.grid = w.mat, psi = psi, z = z, mu = mu, y = y, 
         omega = omega, n.val = n.val, d = d, D = D, Q = Q, X.mu.re = X.mu.re, X.re = X.re, 
         X.omega.re = X.omega.re, X.w = X.w, gamma.star = gamma.star, 
         delta.star = delta.star, beta.star = beta.star)
  )
}
