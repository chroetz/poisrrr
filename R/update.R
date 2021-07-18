svdK <- function (x, K) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  np <- min(n, p)
  u <- matrix(0, n, np)
  vt <- matrix(0, np, p)
  res <- .Internal(La_svd("S", x, double(np), u, vt))
  res$d <- res$d[seq_len(K)]
  res$u <- res$u[, seq_len(K), drop = FALSE]
  res$v <- t(res$vt[seq_len(K), , drop = FALSE])
  res$vt <- NULL
  res
}

svdK_old <- function(x, K) {
  res <- svd(x, K, K)
  res$d <- res$d[1:K]
  res
}


update_f <- function(Y, K, theta, step_size) {
  J <- nrow(Y)
  L <- ncol(Y)
  th_rm <- rowMeans(theta)
  p <- theta - th_rm
  th_cm <- colMeans(p)
  p <- p - rep(th_cm, each=J)
  svd_res <- svdK(p, K)
  b <- svd_res$u
  d <- svd_res$d
  X <- b %*% diag(d, K, K)
  X1 <- cbind(1, X)
  f_new <- cbind(th_cm, svd_res$v)
  for (l in 1:L) {
    f_new[l, ] <- poi_fit1_newton_linesearch(
      X1, Y[, l], offset=th_rm, coef=f_new[l, ], step_size=step_size)
  }
  cbind(X1, th_rm) %*% t(cbind(f_new, 1))
}

update_b <- function(Y, K, theta, step_size) {
  J <- nrow(Y)
  L <- ncol(Y)
  th_rm <- rowMeans(theta)
  p <- theta - th_rm
  th_cm <- colMeans(p)
  p <- p - rep(th_cm, each=J)
  svd_res <- svdK(p, K)
  d <- svd_res$d
  f <- svd_res$v
  X <- f %*% diag(d, K, K)
  X1 <- cbind(1, X)
  b_new <- cbind(th_rm, svd_res$u)
  for (j in 1:J) {
    b_new[j, ] <- poi_fit1_newton_linesearch(
      X1, Y[j, ], offset=th_cm, coef=b_new[j, ], step_size=step_size)
  }
  cbind(1, b_new) %*% t(cbind(th_cm, X1))
}

theta2plist <- function(theta, K) {
  th_m <- mean(theta)
  p <- theta - th_m
  th_rm <- rowMeans(p)
  p <- p - th_rm
  th_cm <- colMeans(p)
  p <- p - rep(th_cm, each=nrow(theta))
  svd_res <- svdK(p, K)
  list(th_m=th_m,
       th_rm=th_rm,
       th_cm=th_cm,
       u=svd_res$u,
       d=svd_res$d,
       v=svd_res$v)
}

plist2theta <- function(plist) {
  K <- length(plist$d)
  J <- length(plist$th_rm)
  theta <-
    plist$u %*% diag(plist$d, K, K) %*% t(plist$v) +
    plist$th_rm + rep(plist$th_cm, each=J) + plist$th_m
  theta
}

loglik <- function(Y, theta) {
  mean(-exp(theta) + Y*theta)
}

loglik_plist_list <- function(plist_list, Y) {
  sapply(plist_list, \(plist) loglik(Y, plist2theta(plist)))
}

