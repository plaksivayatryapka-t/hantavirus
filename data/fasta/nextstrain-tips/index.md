source: https://nextstrain.org/groups/hodcroftlab/andv
fasta_source: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
inputs:
- data/nextstrain/groups/hodcroftlab/andv/S.json
- data/nextstrain/groups/hodcroftlab/andv/M.json
- data/nextstrain/groups/hodcroftlab/andv/L.json
outputs:
- S.fasta
- M.fasta
- L.fasta
- manifest.tsv

FASTA records are fetched from NCBI nuccore accessions found on leaf nodes in the Nextstrain Auspice JSON files.
