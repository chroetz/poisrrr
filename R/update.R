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

update_f <- function(Y, K, theta) {
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
    f_new[l, ] <- step1_newton_pois(
      X1, Y[, l], offset=th_rm, coef=f_new[l, ])
  }
  cbind(X1, th_rm) %*% t(cbind(f_new, 1))
}

update_b <- function(Y, K, theta) {
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
    b_new[j, ] <- step1_newton_pois(
      X1, Y[j, ], offset=th_cm, coef=b_new[j, ])
  }
  cbind(1, b_new) %*% t(cbind(th_cm, X1))
}
