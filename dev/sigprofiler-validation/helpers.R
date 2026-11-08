# Helpers shared by compare.R, compare_tcga.R and compare_pcawg.R: reading
# SigProfilerMatrixGenerator (SPMG) output and comparing it with musicatk.
# Part of the validation for issue #218; not part of the package.

# SPMG ID83 name to musicatk IND83 name: 1:Del:C:5 is DEL_C_1_5+,
# 2:Ins:R:0 is INS_repeats_2_0 and 5:Del:M:5 is DEL_MH_5+_5+
.id83_to_musicatk <- function(x) {
  if (length(x) == 0) return(character(0))
  p <- do.call(rbind, strsplit(x, ":"))
  size <- ifelse(p[, 1] == "5", "5+", p[, 1])
  n <- ifelse(p[, 4] == "5", "5+", p[, 4])
  type <- toupper(p[, 2])
  kind <- c(R = "repeats", M = "MH")[p[, 3]]
  ifelse(p[, 1] == "1", paste(type, p[, 3], "1", n, sep = "_"),
    paste(type, kind, size, n, sep = "_")
  )
}

# An SPMG count matrix (SBS96, DBS78 or ID83) with musicatk motif names, in
# the order of 'motifs'. 'spmg_dir' is SPMG's output_directory. A missing
# matrix is an error unless 'required' is FALSE (SPMG writes no DBS matrix
# when there are no DBS), when NULL is returned with a message.
.read_spmg_matrix <- function(spmg_dir, project, ext, motifs,
                              required = TRUE) {
  path <- file.path(
    spmg_dir, sub("[0-9]+$", "", ext), paste0(project, ".", ext, ".all")
  )
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("SPMG wrote no ", ext, " matrix: ", path)
    message("SPMG wrote no ", ext, " matrix for ", project)
    return(NULL)
  }
  m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  if (ext == "ID83") {
    rownames(m) <- .id83_to_musicatk(rownames(m))
    if (!setequal(rownames(m), motifs)) stop("ID83 names don't match IND83")
    m[motifs, , drop = FALSE]
  } else {
    .match_count_rows(m, motifs, ext)
  }
}

# SPMG's per-variant indel output (seqInfo = True): sample, chrom, position
# of the base before the indel, type, inserted or deleted sequence, and ID83
# bin as a musicatk name. SPMG writes the alleles with an anchor base (VCF)
# or a "-" (MAF) first, so the sequence is the longer allele without it
.read_spmg_indels <- function(spmg_dir) {
  files <- list.files(file.path(spmg_dir, "vcf_files", "ID"),
    full.names = TRUE
  )
  rows <- lapply(files[file.size(files) > 0], read.delim,
    header = FALSE, colClasses = "character", quote = ""
  )
  d <- do.call(rbind, rows)
  if (is.null(d)) {
    return(data.frame(
      sample = character(0), chrom = character(0), pos = integer(0),
      type = character(0), seq = character(0), spmg = character(0)
    ))
  }
  label <- sub("^[A-Z]:", "", d$V4)
  longer <- ifelse(nchar(d$V5) > nchar(d$V6), d$V5, d$V6)
  data.frame(
    sample = d$V1, chrom = d$V2, pos = as.integer(d$V3),
    type = ifelse(grepl(":Ins:", label), "INS", "DEL"),
    seq = substring(longer, 2), spmg = .id83_to_musicatk(label)
  )
}

# Compare two motifs x samples count matrices with the same rows. Samples
# missing from one side count as zero there and are reported.
.compare_counts <- function(a, b, name, a_name = "musicatk", b_name = "spmg") {
  if (is.null(b)) b <- a[, 0, drop = FALSE]
  if (is.null(a)) a <- b[, 0, drop = FALSE]
  samples <- union(colnames(a), colnames(b))
  full <- function(m) {
    out <- matrix(0, nrow(a), length(samples),
      dimnames = list(rownames(a), samples)
    )
    out[, colnames(m)] <- m[rownames(a), , drop = FALSE]
    out
  }
  fa <- full(a)
  fb <- full(b)
  diff <- which(fa != fb, arr.ind = TRUE)
  cells <- data.frame(
    table = rep(name, nrow(diff)),
    sample = samples[diff[, "col"]],
    motif = rownames(fa)[diff[, "row"]],
    a = fa[diff], b = fb[diff]
  )
  names(cells)[4:5] <- c(a_name, b_name)
  summary <- data.frame(
    table = name, samples = length(samples),
    a_total = sum(fa), b_total = sum(fb),
    samples_identical = sum(colSums(fa != fb) == 0),
    samples_only_a = length(setdiff(samples, colnames(b))),
    samples_only_b = length(setdiff(samples, colnames(a))),
    cells_different = nrow(diff),
    abs_difference = sum(abs(fa - fb))
  )
  names(summary) <- sub("_a$", paste0("_", a_name), names(summary))
  names(summary) <- sub("^a_", paste0(a_name, "_"), names(summary))
  names(summary) <- sub("_b$", paste0("_", b_name), names(summary))
  names(summary) <- sub("^b_", paste0(b_name, "_"), names(summary))
  list(summary = summary, cells = cells)
}

# Per-variant IND83 bins of the indels in musica 'm' from musicatk's
# classifier, keyed like .read_spmg_indels(): position of the base before
# a deletion, or of the base an insertion follows, and the inserted or
# deleted sequence
.musicatk_indels <- function(m, g) {
  var <- as.data.frame(variants(m))
  indels <- var[var$Variant_Type %in% c("INS", "DEL"), ]
  data.frame(
    sample = as.character(indels$sample),
    chrom = sub("^chr", "", indels$chr),
    pos = ifelse(indels$Variant_Type == "DEL",
      indels$start - 1L, indels$start
    ),
    type = indels$Variant_Type,
    seq = ifelse(indels$Variant_Type == "INS", indels$alt, indels$ref),
    musicatk = .classify_indels(indels, g)
  )
}

# Join per-variant bins from two sources by sample, chrom, pos, type and
# sequence, so that two indels at one position (a multiallelic site) don't
# pair up crosswise. Records repeated within a source are paired in order
# (first with first) rather than crosswise, and counted as duplicates.
# Returns the joined table, counts by type (INS, DEL), and the disagreeing
# bins as a cross-table of (bin in x, bin in y, n).
.join_indels <- function(x, y, x_col, y_col) {
  keys <- c("sample", "chrom", "pos", "type", "seq")
  occurrence <- function(d) {
    key <- do.call(paste, d[keys])
    stats::ave(seq_along(key), key, FUN = seq_along)
  }
  x$occurrence <- occurrence(x)
  y$occurrence <- occurrence(y)
  j <- merge(x, y, by = c(keys, "occurrence"), all = TRUE)
  in_x <- !is.na(j[[x_col]])
  in_y <- !is.na(j[[y_col]])
  both <- in_x & in_y
  same <- both & j[[x_col]] == j[[y_col]]
  counts <- do.call(rbind, lapply(c("INS", "DEL"), function(t) {
    k <- j$type == t
    data.frame(
      type = t, matched = sum(k & both), agree = sum(k & same),
      only_x = sum(k & in_x & !in_y), only_y = sum(k & in_y & !in_x),
      duplicates_x = sum(k & in_x & j$occurrence > 1),
      duplicates_y = sum(k & in_y & j$occurrence > 1)
    )
  }))
  disagree <- j[both & !same, ]
  tab <- data.frame(x = character(0), y = character(0), n = integer(0))
  if (nrow(disagree) > 0) {
    tab <- as.data.frame(table(disagree[[x_col]], disagree[[y_col]]),
      stringsAsFactors = FALSE
    )
  }
  names(tab) <- c(x_col, y_col, "n")
  list(counts = counts, disagree = tab[tab$n > 0, ], joined = j)
}

# Delete a working folder and fail loudly if anything is left
.remove_dir <- function(dir) {
  unlink(dir, recursive = TRUE)
  if (dir.exists(dir)) stop("Could not delete ", dir, "; delete it by hand")
  invisible(TRUE)
}

# Run expr and keep its messages and warnings without printing them
.quietly <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, paste("Warning:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = trimws(msgs))
}

# Run SPMG on one prepared input folder with run_spmg.py --dir
.run_spmg <- function(python, input_dir, project, genome, out_dir) {
  script <- file.path("dev", "sigprofiler-validation", "run_spmg.py")
  status <- system2(python,
    c(script, "--dir", input_dir, project, genome, out_dir),
    stdout = FALSE, stderr = FALSE
  )
  if (status != 0) stop("SPMG failed on ", input_dir)
  invisible(out_dir)
}
