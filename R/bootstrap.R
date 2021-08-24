#' Estimate the parameters of the INAR model for one token
#'
#' @param y counts of a token for each point in time
#' @param theta parameter vector for the token for each point in time
estim_INAR_token <- function(y, theta) {
  m  <- length(y) # number of points in time
  mu <- exp(theta)
  delta <- mu[2:m] - mu[1:(m-1)]
  b  <- 10^8 # Value is arbitrarily set.
  a <- max(c(0, delta))
  minf <- stats::optimize(inar_objective, lower=a, upper=b, y=y, mu=mu)
  gamma <- minf$minimum
  eta <- mu[2:m] / (mu[1:(m-1)] + gamma)
  return(c(gamma, eta))
}

estim_INAR_group <- function(y, th) {
  m <- ncol(y)
  if (m==1) return(list(mu = exp(th)))

  gamma_eta <- matrix(NA_real_, nrow=ncol(th), ncol=nrow(th))
  dimnames(gamma_eta) <- dimnames(th)[2:1]
  for (i in 1:nrow(y)) {
    gamma_eta[,i] <- estim_INAR_token(y[i,], th[i,])
  }
  eta  <- t(gamma_eta[-1,])
  gamma <- gamma_eta[1,]

  list(
    gamma = gamma,
    eta = eta,
    omega = eta * gamma,
    mu = exp(th))
}

#' Estimate INAR parameters.
#'
#' @param Y A n x d matrix containing token counts for n different tokens and d
#'   different documents.
#' @param theta A n x d matrix. The logs of the Poisson parameters.
#' @param group A character vector of length n for grouping.
#' @param time A numeric vector of length n to order elements by time.
#' @param row_ids A index vector of length n mapping the ordering of group and
#'   time to the ordering in the rows of theta.
#' @return A named list with INAR parameters for each group.
#' @export
estim_INAR <- function(Y, theta,
                       group=rep("default", nrow(Y)),
                       time=seq_len(group),
                       row_ids=rownames(Y)) {
  groups <- unique(group)
  boot_param <- lapply(groups, \(g) {
    sel <- which(group == g)
    tm <- time[sel]
    idx <- sel[order(tm)]
    y_group <- Y[, row_ids[idx], drop=FALSE]
    th_group <- theta[, row_ids[idx], drop=FALSE]
    estim_INAR_group(y_group, th_group)
  })
  names(boot_param) <- groups
  boot_param
}

#' Create bootstrap samples form INAR parameters.
#'
#' @param Y A n x d matrix containing token counts for n different tokens and d
#'   different documents.
#' @param boot_param .
#' @return .
#' @export
boot_sample <- function(Y, boot_param) {
  n <- nrow(Y)
  Y_new <- matrix(NA_real_, nrow=n, ncol=ncol(Y))
  dimnames(Y_new) <- dimnames(Y)

  for (p in boot_param) {
    ids <- colnames(p$mu)
    m <- length(ids)

    Y_new[, ids[1]] <- stats::rpois(n, p$mu[, ids[1]])
    if (m==1) next

    pois_omega <- p$omega
    pois_omega[] <- stats::rpois(length(p$omega), p$omega)
    for (t in 2:m) {
      id_jt <- ids[t]
      xi <- stats::rbinom(nrow(Y_new), size=Y_new[, ids[t-1]], prob=p$eta[, id_jt])
      z_ijt <- pois_omega[, id_jt]
      Y_new[, id_jt] <- xi + z_ijt
    }
  }
  Y_new
}

