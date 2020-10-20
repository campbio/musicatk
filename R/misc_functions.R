.kld <- function(a, b) {
  return(sum(a * log2(a / b)))
}

#' Calculates 1 - Jensen-Shannon Divergences between all pairs of columns 
#' between two matrices
#'
#' @param p First matrix
#' @param q Second matrix
#' @param epsilon Number to add to all probabilities. Default \code{0.0000001}.
#' @return Returns matrix of 1 - Jensen-Shannon Divergences
#' @keywords internal
.jsd <- function(p, q, epsilon = 0.0000001) {
  # Add small value to handle zeros and then renormalize using prop.table
  p <- prop.table(p + epsilon, margin = 2)
  q <- prop.table(q + epsilon, margin = 2)

  res <- matrix(nrow = ncol(p), ncol = ncol(q), dimnames = list(colnames(p),
                                                                colnames(q)))
  for (i in seq_len(ncol(p))) {
    for (j in seq_len(ncol(q))) {
      m <- (p[,i] + q[,j]) / 2
      res[i, j] <- 1 - (0.5 * .kld(p[,i], m) + 0.5 * .kld(q[,j], m))
    }
  }
  return(res)
}


# Calculates distance based off of cosine similarity
.cosineDist <- function(x) {
  y <- (1 - .cosine(x, x)) / 2
  return(stats::as.dist(y))
}

# Calculates cosine similarities between 2 matrices where columns are sigs
# and rows are mutation motifs.
.cosine <- function(x, y) {
  nX <- ncol(x)
  nY <- ncol(y)
  if(nrow(x) != nrow(y)) {
    stop("The number of rows in 'x' and 'y' must be the same.")
  }
  temp <- t(cbind(x, y))
  res <- temp %*% t(temp) / (sqrt(rowSums(temp^2) %*% t(rowSums(temp^2))))
  return(res[1:nX, -(1:nX)])
}


