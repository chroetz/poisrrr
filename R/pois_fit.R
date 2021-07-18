# make one step of newton method for fitting to a Poisson model
step1_newton_pois <- function(X, y, offset=0, coef) {
  eta <- X %*% coef + offset # eta < -708.3964 or eta  > 709.7827 makes mu Inf
  dim(eta) <- NULL
  mu <- exp(eta)
  grad <- (y-mu) %*% X
  dim(grad) <- NULL
  hess <- crossprod(X*mu, X)
  hess_stable <- hess + diag(nrow(hess))*max(abs(hess))*1e-14
  s <- .Internal(La_solve(hess_stable, grad, .Machine$double.eps))

  stp_sz <- 1 # start with step size 1
  ll_old <- mean(-exp(eta) + y*eta)
  while (TRUE) {
    res <- coef + stp_sz * s
    eta_new <- X %*% res + offset
    dim(eta_new) <- NULL
    ll_new <- mean(-exp(eta_new) + y*eta_new)
    if (ll_new >= ll_old) break
    if (stp_sz < 1e-3) {
      warning("Could not find improvment. Returning old coefficient.")
      return(coef)
    }
    stp_sz <- stp_sz / 2
  }

  res
}
