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

#' Rotate a point configuration to fit another one.
#'
#' @param target A numeric n x d matrix.
#' @param x A numeric n x d matrix.
#' @return A n x d matrix containing the rotated points of x which minimize the
#'   squared Euclidean distance to target.
#'
#' @export
rotate_fit <- function(x, target) {
  n1 <- ncol(target)-1
  # optimize wrt det 1 matrices
  gmin_pos <- stats::optim(
    rep(pi, (n1+1)*n1/2),
    rotate_fit_objective,
    method = "L-BFGS-B",
    lower=0, upper=2*pi,
    v=target, v_b=x)
  # optimize wrt det -1 matrices
  v_b_neg <- x %*% diag(c(-1, rep(1, n1)))
  gmin_neg <- stats::optim(
    rep(pi, (n1+1)*n1/2),
    rotate_fit_objective,
    method = "L-BFGS-B",
    lower=0, upper=2*pi,
    v=target, v_b=v_b_neg)
  if (gmin_pos$value <= gmin_neg$value) {
    return(rotate(gmin_pos$par, x))
  }
  return(rotate(gmin_neg$par, v_b_neg))
}

