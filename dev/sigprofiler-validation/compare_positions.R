# Compare musicatk IND83 with SigProfilerMatrixGenerator (SPMG) ID83 on
# indels written at different positions along a repeat or a microhomology
# (issue #235). Part of the SPMG validation; not part of the package.
#
# Usage, from the package root:
#   Rscript dev/sigprofiler-validation/compare_positions.R <out_dir> <python>
#
# Generates indels on hg38 chr1:10-30 Mb: random 1-8 bp deletions, random
# 1-8 bp tandem duplications, and one unit inserted into or deleted from
# 2-5 bp microsatellites. Each event is written as three VCF samples: left-
# aligned ("left"), at a random position in between ("mid") and right-
# aligned ("right"), and counted by both tools. Writes bins.csv to <out_dir>.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compare_positions.R <out_dir> <python>")
}
out_dir <- args[1]
python <- args[2]
suppressMessages(devtools::load_all(".", quiet = TRUE))
source(file.path("dev", "sigprofiler-validation", "helpers.R"))
g <- select_genome("hg38")
set.seed(235)

chr <- "chr1"
region_start <- 10000000L
region_end <- 30000000L
ref_seq <- as.character(BSgenome::getSeq(g, chr, region_start, region_end))
bases_at <- function(p, q = p) {
  substring(ref_seq, p - region_start + 1L, q - region_start + 1L)
}

# Events are in VCF anchor form: an insertion adds 'unit' after base p, a
# deletion removes bases p + 1 to p + nchar(unit)
shift_left <- function(type, p, unit) {
  len <- nchar(unit)
  repeat {
    if (type == "INS") {
      if (bases_at(p) != substr(unit, len, len)) break
      unit <- paste0(bases_at(p), substr(unit, 1, len - 1))
    } else {
      if (bases_at(p) != bases_at(p + len)) break
      unit <- bases_at(p, p + len - 1)
    }
    p <- p - 1L
  }
  list(p = p, unit = unit)
}
# Every position of the event, from the left-aligned one rightwards
positions <- function(type, p, unit) {
  len <- nchar(unit)
  out <- list(list(p = p, unit = unit))
  repeat {
    if (type == "INS") {
      if (bases_at(p + 1) != substr(unit, 1, 1)) break
      unit <- paste0(substr(unit, 2, len), bases_at(p + 1))
    } else {
      if (bases_at(p + 1) != bases_at(p + len + 1)) break
      unit <- bases_at(p + 2, p + len + 1)
    }
    p <- p + 1L
    out[[length(out) + 1]] <- list(p = p, unit = unit)
  }
  out
}

# Random positions spaced at least 500 bp apart, with no N nearby
n_each <- 1000
random_positions <- function(n) {
  pos <- sort(sample(
    seq(region_start + 1000L, region_end - 1000L, by = 500L), n
  ))
  pos[!grepl("N", bases_at(pos - 300L, pos + 300L))]
}
random_len <- function(n) sample(1:8, n, replace = TRUE)

del_pos <- random_positions(n_each)
del_len <- random_len(length(del_pos))
dup_pos <- random_positions(n_each)
dup_len <- random_len(length(dup_pos))

# Microsatellites of 2-5 bp units (not homopolymers) with three or more
# copies, in the first 5 Mb of the region
chunk <- substr(ref_seq, 1, 5e6)
hits <- gregexpr("([ACGT]{2,5})\\1{2,}", chunk, perl = TRUE)[[1]]
cap_start <- attr(hits, "capture.start")[, 1]
cap_len <- attr(hits, "capture.length")[, 1]
str_unit <- substring(chunk, cap_start, cap_start + cap_len - 1)
usable <- which(hits > 1000 & !grepl("^(.)\\1+$", str_unit, perl = TRUE))
usable <- sample(usable, min(n_each, length(usable)))
str_offset <- vapply(usable, function(i) {
  sample(0:(attr(hits, "match.length")[i] - cap_len[i]), 1)
}, integer(1))
str_pos <- region_start + as.integer(hits[usable]) - 2L + str_offset
str_ins <- runif(length(usable)) < 0.5

events <- data.frame(
  source = c(
    rep("random_del", length(del_pos)), rep("random_dup", length(dup_pos)),
    ifelse(str_ins, "str_ins", "str_del")
  ),
  type = c(
    rep("DEL", length(del_pos)), rep("INS", length(dup_pos)),
    ifelse(str_ins, "INS", "DEL")
  ),
  p = c(del_pos, dup_pos, str_pos),
  len = c(del_len, dup_len, cap_len[usable])
)
events$unit <- bases_at(events$p + 1L, events$p + events$len)
events <- events[!grepl("N", events$unit), ]

rows <- lapply(seq_len(nrow(events)), function(i) {
  e <- events[i, ]
  left <- shift_left(e$type, e$p, e$unit)
  all <- positions(e$type, left$p, left$unit)
  n <- length(all)
  mid <- if (n >= 3) all[[sample(2:(n - 1), 1)]] else all[[1]]
  picks <- list(left = all[[1]], mid = mid, right = all[[n]])
  do.call(rbind, lapply(names(picks), function(s) {
    r <- picks[[s]]
    anchor <- bases_at(r$p)
    with_unit <- paste0(anchor, r$unit)
    data.frame(
      event = i, source = e$source, type = e$type, len = e$len,
      n_positions = n, first = all[[1]]$p, last = all[[n]]$p,
      sample = s, chr = chr, pos = r$p,
      ref = if (e$type == "INS") anchor else with_unit,
      alt = if (e$type == "INS") with_unit else anchor
    )
  }))
})
ev <- do.call(rbind, rows)

# Keep events whose positions are more than 50 bp from any other event's,
# so each SPMG record maps to one event
spans <- unique(ev[, c("event", "first", "last")])
spans <- spans[order(spans$first), ]
gap_before <- c(Inf, spans$first[-1] - spans$last[-nrow(spans)])
gap_after <- c(gap_before[-1], Inf)
ev <- ev[ev$event %in% spans$event[gap_before > 50 & gap_after > 50], ]

input_dir <- file.path(out_dir, "input")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
header <- c("##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)
for (s in c("left", "mid", "right")) {
  d <- ev[ev$sample == s, ]
  d <- d[order(d$pos), ]
  records <- paste(d$chr, d$pos, ".", d$ref, d$alt, ".", "PASS", ".",
    sep = "\t"
  )
  writeLines(c(header, records), file.path(input_dir, paste0(s, ".vcf")))
}

# musicatk bins, from the same records
dt <- data.frame(
  chr = ev$chr, start = ev$pos, end = ev$pos + nchar(ev$ref) - 1L,
  ref = ev$ref, alt = ev$alt, sample = ev$sample
)
m <- suppressMessages(create_musica_from_variants(dt, g, verbose = FALSE))
mus <- .musicatk_indels(m, g)
ev$musicatk <- mus$musicatk[match(
  paste(ev$sample, ev$pos), paste(mus$sample, mus$pos)
)]

# SPMG bins, matched to events by sample and span of positions
spmg_dir <- file.path(out_dir, "spmg")
.run_spmg(python, input_dir, "positions", "GRCh38", spmg_dir)
sp <- .read_spmg_indels(spmg_dir)
spans <- unique(ev[, c("event", "first", "last")])
spans <- spans[order(spans$first), ]
k <- findInterval(sp$pos, spans$first)
k[k == 0] <- NA
k[!is.na(k) & sp$pos > spans$last[pmax(k, 1)]] <- NA
sp$event <- spans$event[k]
ev$spmg <- sp$spmg[match(
  paste(ev$sample, ev$event), paste(sp$sample, sp$event)
)]
write.csv(ev, file.path(out_dir, "bins.csv"), row.names = FALSE)

shifted <- ev$n_positions >= 3
cat("Events:", length(unique(ev$event)), "\n")
cat("Events with a middle position:", length(unique(ev$event[shifted])),
  "\n"
)
cat("SPMG records not matched to an event:", sum(is.na(sp$event)), "\n")
cat("\nmusicatk agrees with SPMG, by sample:\n")
print(tapply(ev$musicatk == ev$spmg, ev$sample, sum))
cat("\nmusicatk (any position) agrees with SPMG on the left-aligned record:\n")
left <- ev[ev$sample == "left", ]
print(tapply(
  ev$musicatk == left$spmg[match(ev$event, left$event)], ev$sample, sum
))
cat("\nEvents whose middle record SPMG bins differently from the left:\n")
mid <- ev[ev$sample == "mid" & shifted, ]
differs <- mid$spmg != left$spmg[match(mid$event, left$event)]
print(table(mid$source[differs]))
