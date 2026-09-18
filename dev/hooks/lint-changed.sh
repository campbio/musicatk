#!/usr/bin/env bash
# PostToolUse hook: report-only lint of a just-edited R file.
#
# Reports lints; never rewrites the file. Auto-styling on edit buries the real
# change in formatting noise whenever a lint backlog exists — run `make style`
# deliberately, as its own commit, instead.
#
# Reads the Claude Code hook payload on stdin and exits 0 unconditionally: a
# lint finding is information, not a reason to fail the edit.

set -uo pipefail

payload="$(cat)"

# Extract .tool_input.file_path without hard-depending on jq.
if command -v jq >/dev/null 2>&1; then
  file="$(printf '%s' "$payload" | jq -r '.tool_input.file_path // empty')"
elif command -v python3 >/dev/null 2>&1; then
  file="$(printf '%s' "$payload" | python3 -c \
    'import json,sys
try:
    print(json.load(sys.stdin).get("tool_input", {}).get("file_path", ""))
except Exception:
    print("")')"
else
  exit 0
fi

[ -n "$file" ] || exit 0
[ -f "$file" ] || exit 0

case "$file" in
  *.R|*.r) ;;
  *) exit 0 ;;
esac

command -v Rscript >/dev/null 2>&1 || exit 0

Rscript -e '
  args <- commandArgs(trailingOnly = TRUE)
  f <- args[1]
  if (!requireNamespace("lintr", quietly = TRUE)) quit(status = 0)
  l <- tryCatch(lintr::lint(f), error = function(e) NULL)
  if (is.null(l) || length(l) == 0) quit(status = 0)
  cat(sprintf("lintr: %d lint(s) in %s\n", length(l), basename(f)))
  print(l)
' "$file" 2>/dev/null

exit 0
