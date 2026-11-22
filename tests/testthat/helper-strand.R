# musica with "chr" stripped from the chromosome names of the given samples
strip_chr <- function(musica, samples = NULL) {
  dat <- variants(musica)
  rows <- if (is.null(samples)) TRUE else dat$sample %in% samples
  dat$chr[rows] <- sub("^chr", "", dat$chr[rows])
  variants(musica) <- dat
  musica
}
