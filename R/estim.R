estim_essential <- function(Y, K, iter=3, log_zero=-1.5,
                            step_size=1) {
  theta <- log(Y)
  theta[theta == -Inf] <- log_zero
  r <- theta2plist(theta, K)
  th <- plist2theta(r)
  for (i in seq_len(iter)) {
    th <- update_f(Y, K, th, step_size)
    th <- update_b(Y, K, th, step_size)
  }
  th
}

estim <- function(
    Y, K, iter=3,
    log_zero=-1.5, step_size=1) {
  pt <- proc.time()
  plist_list <- list()
  theta <- log(Y)
  theta[!is.finite(theta)] <- log_zero
  r <- theta2plist(theta, K)
  th <- plist2theta(r)
  plist_list[[length(plist_list)+1]] <- theta2plist(th, K)
  cat("0%")
  elapsed <- double(2*iter+1)
  elapsed[1] <- (proc.time() - pt)[3]
  for (i in seq_len(iter)) {
    elapsed[i*2] <- system.time({
      th <- update_f(Y, K, th, step_size)
    })[3]
    plist_list[[length(plist_list)+1]] <- theta2plist(th, K)
    elapsed[i*2+1] <- system.time({
      th <- update_b(Y, K, th, step_size)
    })[3]
    plist_list[[length(plist_list)+1]] <- theta2plist(th, K)
    cat("\r", round(i/iter*100), "%", sep="")
  }
  cat(" Done.\n")
  kind <- c("init", rep(c("update_f", "update_b"), times=iter))
  tibble(
    plist = plist_list,
    loop_num = c(0, rep(1:iter, each=2)),
    kind = kind,
    step = 0:(2*iter),
    elapsed = elapsed)
}

estim_tb <- function(opts, Y) {
  results <- list()
  for (i in seq_len(nrow(opts))) {
    cat(i, "/", nrow(opts), ":\n", sep="")
    res <- estim(Y, opts$K[[i]],
                 log_zero = opts$log_zero[[i]],
                 iter=opts$iter[[i]],
                 step_size=opts$step_size[[i]])
    results[[length(results)+1]] <-
      res %>% mutate(loglik = loglik_plist_list(plist, Y))
  }
  opts$results <- results
  opts %>% rowid_to_column(var = "opt_id")
}

eval_estim <- function(res, Y) {
  res %>%
    mutate(loglik = loglik_plist_list(plist, Y)) ->
    iterations
  row_best <- which.max(iterations$loglik)
  list(
    iterations = iterations,
    elapsed = sum(iterations$elapsed),
    loglik = iterations$loglik[[row_best]],
    plist = iterations$plist[[row_best]])
}
