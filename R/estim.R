estim <- function(Y, K, tol=1e-5, verbose=TRUE) {
  theta <- log(Y)
  theta[!is.finite(theta)] <- -2 # set log(0) to -2
  r <- theta2plist(theta, K)
  th <- plist2theta(r)
  ll_old <- loglik(Y, th)
  i <- 1
  while(TRUE) {
    th <- update_f(Y, K, th)
    th <- update_b(Y, K, th)
    ll_new <- loglik(Y, th)
    rel_improve <- (ll_new - ll_old) / abs(ll_old)
    if (verbose) cat(sprintf("step %d: relative improvement %f\n", i, rel_improve))
    if (rel_improve < tol) break
    ll_old <- ll_new
    i <- i+1
  }
  th
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

