rot_mat_hyperplane <- function(alpha, i, j, n) {
  mat <- diag(n)
  mat[i,i] <- cos(alpha)
  mat[j,j] <- cos(alpha)
  mat[i,j] <- -sin(alpha)
  mat[j,i] <- sin(alpha)
  mat
}

rotate <- function(alpha, x) {
  n <- ncol(x)
  z <- 0
  for (i in 1:(n-1)) for (j in (i+1):n) {
    z <- z + 1
    x <- x %*% rot_mat_hyperplane(alpha[z], i, j, n)
  }
  x
}

rotate_fit_objective <- function(alpha, v, v_b) {
  sum((v - rotate(alpha, v_b))^2)
}

#' @export
rotate_fit <- function(v_b, v) {
  n1 <- ncol(v)-1
  gmin <- optim(rep(pi, (n1+1)*n1/2), rotate_fit_objective,
                method = "L-BFGS-B", lower=0, upper=2*pi, v=v, v_b=v_b)
  rotate(gmin$par , v_b)
}

