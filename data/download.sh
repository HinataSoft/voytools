#!/usr/bin/env bash
# Re-download the transliteration files from voynich.nu into this folder.
# Needs network access to www.voynich.nu (in a cloud session: environment network
# access "Custom" with www.voynich.nu allowed, or "Full").
set -euo pipefail
cd "$(dirname "$0")"
base="https://www.voynich.nu/data"
for f in ZL3b-n.txt GC2a-n.txt RF1b-er.txt IT2a-n.txt; do
  echo "GET $f"; curl -fsSL -o "$f" "$base/$f"
done
echo "GET LSI_ivtff_0d.txt"; curl -fsSL -o LSI_ivtff_0d.txt "$base/beta/LSI_ivtff_0d.txt"
for f in ZL3b-n.txt GC2a-n.txt RF1b-er.txt IT2a-n.txt LSI_ivtff_0d.txt; do
  head -c 7 "$f" | grep -q '#=IVTFF' || { echo "ERROR: $f is not an IVTFF file" >&2; exit 1; }
done
echo "OK"
