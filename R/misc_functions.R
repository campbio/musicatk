kld <- function(a, b) {
  return(sum(a * log2(a / b)))
}

#' Compare two vectors similarity based on Jensen-Shannon Divergence
#'
#' @param p First vector
#' @param q Second vector
#' @return Returns Jensen-Shannon Divergence of the two vectors
#' @examples
#' p <- c(0.2, 0.3, 0.4, 0.5, 0.6)
#' q <- c(0.6, 0.5, 0.4, 0.3, 0.2)
#' jsd(p, q)
#' @export
jsd <- function(p, q) {
  epsilon <- 0.0000001
  p <- p + epsilon
  q <- q + epsilon
  m <- (p + q) / 2
  jsd <- 0.5 * kld(p, m) + 0.5 * kld(q, m)
  return(jsd)
}

