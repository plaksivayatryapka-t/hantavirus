# Data Sources and Attribution

This repository contains sequence data and derived analysis artifacts for local research use.

## NCBI GenBank / nuccore

FASTA records under `data/fasta/` were retrieved from NCBI GenBank/nuccore accessions.

NCBI states that it places no restrictions on use or distribution of GenBank data, while noting that original submitters may claim patent, copyright, or other intellectual property rights in submitted data and that NCBI cannot assess those claims.

Sources:
- https://www.ncbi.nlm.nih.gov/genbank/about/
- https://www.ncbi.nlm.nih.gov/home/about/policies/

## Nextstrain

Nextstrain/Auspice JSON files under `data/nextstrain/` were retrieved from:

- https://nextstrain.org/groups/hodcroftlab/andv

Derived alignment and visualization artifacts may use these JSON files and NCBI FASTA records as inputs. When sharing figures or analysis results derived from Nextstrain views, credit Nextstrain and the underlying genomic data contributors.

Nextstrain notes that screenshots are licensed under CC-BY and should credit the authors behind the genomic data as well as Nextstrain.

Sources:
- https://docs.nextstrain.org/en/latest/guides/share/download-data.html
- https://nextstrain.org/

## Local Derived Artifacts

Files under `data/alignments/` are locally generated analysis outputs. Each analysis directory should include an `index.md` describing its inputs, method, and outputs.
