#!/usr/bin/env bash
# Run from the repository root BEFORE the first commit and before每次 release.
# Fails loudly if anything identifying or machine-specific is still present.
set -u
PATTERNS=(
  "D:/download"                 # Windows working directory
  "/n/holylfs"                  # cluster path
  "christiani_lab"              # lab directory
  "holy2a24315"                 # hostname
  "rc.fas.harvard.edu"          # cluster domain
  "zhanxiali"                   # user account
  "1999P004935"                 # IRB protocol number
  "setwd("                      # any hard-coded working directory
)
fail=0
for p in "${PATTERNS[@]}"; do
  if grep -rInF -- "$p" . --include="*.R" --include="*.cpp" --include="*.md" --include="*.txt" 2>/dev/null | grep -v scrub_check.sh; then
    echo ">>> FOUND: $p"; fail=1
  fi
done
if grep -rInE '"[A-Za-z]:[/\\]|/Users/|/home/[a-z]+/' . --include="*.R" --include="*.cpp" 2>/dev/null | grep -v scrub_check.sh; then
  echo ">>> FOUND: absolute path literal"; fail=1
fi
[ $fail -eq 0 ] && echo "Clean - no identifying paths or identifiers found." || echo "SCRUB REQUIRED - fix the hits above before committing."
