source: https://nextstrain.org/groups/hodcroftlab/andv/M
method: MAFFT v7.505
segment: M
reason: M encodes the Andes virus glycoprotein precursor, relevant to entry/tropism comparisons.

input:
- Swiss_2026_M from data/fasta/ANDV-Switzerland-Hu-3337-2026/ANDV-Switzerland-Hu-3337-2026.fasta
- NRC-4_2018_M_MN258192.1 from data/fasta/nextstrain-tips/M.fasta
- NRC-3_1997_M_MN258191.1 from data/fasta/nextstrain-tips/M.fasta

output:
- input.fasta
- aligned.fasta
- M_closest_variable_sites.png
- protein_input.fasta
- protein_aligned.fasta
- M_glycoprotein_closest_variable_sites.png
- summary.tsv
- variable_sites.tsv
inferred_m_glycoprotein_orfs:
- Swiss_2026_M: nt 44-3460, 1138 aa
- NRC-4_2018_M_MN258192.1: nt 50-3466, 1138 aa
- NRC-3_1997_M_MN258191.1: nt 48-3464, 1138 aa
