#' Estimate Poisson parameter matrix with reduced rank.
#'
#' @param Y A n x d matrix containing words counts for n different words and d
#'   different documents.
#' @param K A positive integer. The rank of the normalized parameter matrix.
#' @param tol A positive double. Stop iterative optimization of log-likelihood
#'   (ll) if improvement is less than tol*ll.
#' @param max_iter A positive integer. Stop iterative optimization of
#'   log-likelihood after max_iter iterations.
#' @param verbose A logical value. Should progress messages be printed on the
#'   console?
#' @param return_list A logical value or a numeric vector. Should a list of
#'   parameter matrices for each iteration (or all iterations with number in
#'   return_list) be returned or just the last?
#' @return A n x d matrix or a list of such matrices. A matrix contains the logs
#'   of the Poisson parameters.
#' @export
estim <- function(Y, K, verbose=TRUE, tol=1e-5, max_iter=100, return_list=FALSE) {
  theta <- log(Y)
  theta[!is.finite(theta)] <- -2 # set log(0) to -2
  r <- theta2plist(theta, K)
  th <- plist2theta(r)
  if (!isFALSE(return_list)) th_list <- list(th)
  ll_old <- loglik(Y, th)
  i <- 1
  while(TRUE) {
    th <- update_f(Y, K, th)
    th <- update_b(Y, K, th)
    if (isTRUE(return_list) || i %in% return_list) th_list <- c(th_list, list(th))
    ll_new <- loglik(Y, th)
    rel_improve <- (ll_new - ll_old) / abs(ll_old)
    if (verbose)
      cat(paste0(sprintf("#%03d: ",i),
                 "loglik ", formatC(ll_new, format="e", digits=4), ", ",
                 "rel_gain ", formatC(rel_improve, format="e", digits=2),
                 "\n"))
    if (rel_improve < tol) break
    ll_old <- ll_new
    i <- i+1
    if (i > max_iter) break
  }
  if (!isFALSE(return_list)) {
    return(
      lapply(th_list, function(th) {dimnames(th) <- dimnames(Y); th}))
  }
  dimnames(th) <- dimnames(Y)
  th
}

#' Extract parameter vectors from parameter matrix.
#'
#' @param theta A n x d matrix. The logs of the Poisson parameteres.
#' @param K A positive integer. The rank of the normalized parameter matrix.
#' @return A list containing the normalizing values and the SVD of the
#'   normalized matrix.
#' @export
theta2plist <- function(theta, K) {
  th_m <- mean(theta)
  p <- theta - th_m
  th_rm <- rowMeans(p)
  names(th_rm) <- rownames(theta)
  p <- p - th_rm
  th_cm <- colMeans(p)
  names(th_cm) <- colnames(theta)
  p <- p - rep(th_cm, each=nrow(theta))
  svd_res <- svdK(p, K)
  u <- svd_res$u
  rownames(u) <- rownames(theta)
  v <- svd_res$v
  rownames(v) <- colnames(theta)
  list(th_m=th_m,
       th_rm=th_rm,
       th_cm=th_cm,
       u=u,
       d=svd_res$d,
       v=v)
}

#' Combine a list of parameter vectors to the parameter matrix.
#'
#' @param plist A list of parameter vectors as returned by theta2plist().
#' @return A n x d matrix. Contains the logs of the Poisson parameters.
#' @export
plist2theta <- function(plist) {
  K <- length(plist$d)
  J <- length(plist$th_rm)
  theta <-
    plist$u %*% diag(plist$d, K, K) %*% t(plist$v) +
    plist$th_rm + rep(plist$th_cm, each=J) + plist$th_m
  theta
}

#' Log-likelihood
#'
#' @param Y A n x d matrix containing words counts for n different words and d
#'   different documents.
#' @param theta A n x d matrix. The logs of the Poisson parameteres.
#' @return A number. The Log-likelihood.
#' @export
loglik <- function(Y, theta) {
  mean(-exp(theta) + Y*theta)
}

