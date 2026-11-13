# 0004. Store benchmark summaries as data.frames

- **Status:** Proposed
- **Date:** 2026-11-13
- **Deciders:** Joshua D. Campbell (maintainer)

## Context

`full_benchmark` has two summary slots, `method_view_summary` and
`sig_view_summary`, and `single_benchmark` has `single_summary`. All three
are declared as `matrix`. Each summary mixes numbers (counts, reconstruction
error, cosine similarity) with one text row ("Sigs Found" or "Methods
Found"), so the matrices are character, and every number comes back as a
string (#201). Users have to convert with `as.numeric()` before comparing
runs, which is the purpose of the summaries. That matters more once several
discovery methods, including the planned Stan models (#206), are
benchmarked side by side.

Changing a slot's class is an object-model change, which needs an ADR. PR I
of #208 rewrites the summary code for #191 and #198 and rebuilds
`full_benchmark_example`, so the change costs little if it is made now.

## Decision

We will declare `method_view_summary`, `sig_view_summary` and
`single_summary` as `data.frame`.

- `method_view_summary` has one row per benchmark run, with the method ID
  as the row name. `sig_view_summary` has one row per true signature.
  `single_summary` is the run's one-row data.frame.
- Column names are the current row labels ("Num Found", "Mean CS", ...).
  Counts are integer, reconstruction error and cosine similarity are
  numeric, and the signature or method lists are character.
- The replacement methods take a `data.frame` instead of a `matrix`.
- `benchmark()` stops with a message to recreate the object with
  `create_benchmark()` when it is given a `full_benchmark` that still holds
  matrix summaries.

## Consequences

- Summaries can be sorted, filtered and plotted directly, and runs
  compared by column.
- The layout is transposed: runs and signatures are rows, not columns.
  Code that indexed the old matrices (`summary["Mean CS", "nmf8"]`) has to
  change (`summary["nmf8", "Mean CS"]`), and printed output looks
  different. The benchmarking article is updated.
- `full_benchmark` objects saved by older versions can still be loaded and
  read, but can't be extended with `benchmark()`; they have to be rebuilt.
  `full_benchmark_example` is rebuilt by
  `data-raw/make_full_benchmark_example.R`.
- Code that set a summary to a matrix with `method_view_summary<-`,
  `sig_view_summary<-` or `single_summary<-` gets a dispatch error. These
  setters are exported, but nothing in the package or its documentation
  sets them by hand.

## Alternatives considered

- **Keep character matrices** and fix only the warnings and dead checks.
  The numbers would stay strings, which is the problem in #201.
- **Numeric matrices without the text rows.** No class change, but the
  summary would lose the list of matched signatures and methods. Those
  would have to be pulled from each run's comparison table.
- **data.frames in the current orientation** (metrics as rows). Each column
  would still mix numbers and text, so it would be character again.
