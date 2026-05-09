#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out_dir="$repo_root/data/alignments/ANDV-Switzerland-Hu-3337-2026-M-closest"

nt_input="$out_dir/input.fasta"
nt_aligned="$out_dir/aligned.fasta"
protein_input="$out_dir/protein_input.fasta"
protein_aligned="$out_dir/protein_aligned.fasta"

mafft --auto "$nt_input" > "$nt_aligned"
mafft --auto "$protein_input" > "$protein_aligned"

echo "wrote $nt_aligned"
echo "wrote $protein_aligned"
