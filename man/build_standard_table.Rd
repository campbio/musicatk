% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_tables.R
\name{build_standard_table}
\alias{build_standard_table}
\title{Builds count tables using various mutation type schemas}
\usage{
build_standard_table(
  musica,
  g,
  table_name,
  strand_type = NULL,
  overwrite = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{musica}{A \code{\linkS4class{musica}} object.}

\item{g}{A \linkS4class{BSgenome} object indicating which genome
reference the variants and their coordinates were derived from.}

\item{table_name}{Name of standard table to build. One of \code{"SBS96"},
\code{"SBS192"}, \code{"DBS"}, or \code{"Indel"}.}

\item{strand_type}{Strand type to use in SBS192 schema. One of
\code{"Transcript_Strand"} or \code{"Replication_Strand"}.
Only used if \code{table_name = SBS192}.}

\item{overwrite}{If \code{TRUE}, any existing count table with the same
name will be overwritten. If \code{FALSE}, then an error will be thrown
if a table with the same name exists within the \code{musica} object.}

\item{verbose}{Show progress bar for processed samples}
}
\value{
No object will be returned. The count tables will be automatically
added to the \code{musica} object.
}
\description{
Generates count tables for different mutation type schemas which can be 
used as input to the mutational signature discovery or prediction functions.
\code{"SBS96"} generates a table for single base substitutions following the
standard 96 mutation types derived from the trinucleotide context.
\code{"SBS192"} is the 96 mutation type schema with the addition of
transcriptional strand or replication strand information added to each base.
\code{"DBS"} generates a table for the double base substitution schema 
used in COSMIC V3. \code{"Indel"} generates a table for insertions and 
deletions following the schema used in COSMIC V3.
}
\examples{
g <- select_genome("19")

data(musica)
build_standard_table(musica, g, "SBS96", overwrite = TRUE)

data(musica)
annotate_transcript_strand(musica, "19")
build_standard_table(musica, g, "SBS192", "Transcript_Strand")

data(musica)
data(rep_range)
annotate_replication_strand(musica, rep_range)
build_standard_table(musica, g, "SBS192", "Replication_Strand")

data(dbs_musica)
build_standard_table(dbs_musica, g, "DBS", overwrite = TRUE)

data(indel_musica)
build_standard_table(indel_musica, g, table_name = "INDEL")
}
