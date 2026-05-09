#!/usr/bin/env python3
"""Reproduce M-segment nucleotide and protein comparisons for Swiss 2026."""

from __future__ import annotations

import subprocess
from itertools import combinations
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO_ROOT = Path(__file__).resolve().parents[1]
SWISS_FASTA = REPO_ROOT / "data/fasta/ANDV-Switzerland-Hu-3337-2026/ANDV-Switzerland-Hu-3337-2026.fasta"
NEXTSTRAIN_M_FASTA = REPO_ROOT / "data/fasta/nextstrain-tips/M.fasta"
OUT_DIR = REPO_ROOT / "data/alignments/ANDV-Switzerland-Hu-3337-2026-M-closest"

SEQUENCES = {
    "Swiss_2026_M": {
        "source": SWISS_FASTA,
        "match": lambda header: header == "ANDV/Switzerland/Hu-3337/2026_M",
        "label": "Swiss 2026",
    },
    "NRC-4_2018_M_MN258192.1": {
        "source": NEXTSTRAIN_M_FASTA,
        "match": lambda header: header.startswith("MN258192.1 "),
        "label": "NRC-4 2018",
    },
    "NRC-3_1997_M_MN258191.1": {
        "source": NEXTSTRAIN_M_FASTA,
        "match": lambda header: header.startswith("MN258191.1 "),
        "label": "NRC-3 1997",
    },
}

# Nextstrain annotates ANDVsMgp1 on its reference coordinates. The downloaded
# FASTA records have slightly different terminal lengths, so the script infers
# the long M glycoprotein ORF per sequence instead of raw-slicing those coords.
MIN_M_GLYCOPROTEIN_ORF_NT = 3000

NT_INPUT_FASTA = OUT_DIR / "input.fasta"
NT_ALIGNED_FASTA = OUT_DIR / "aligned.fasta"
NT_VARIABLE_PNG = OUT_DIR / "M_closest_variable_sites.png"
AA_INPUT_FASTA = OUT_DIR / "protein_input.fasta"
AA_ALIGNED_FASTA = OUT_DIR / "protein_aligned.fasta"
AA_VARIABLE_PNG = OUT_DIR / "M_glycoprotein_closest_variable_sites.png"
SUMMARY_TSV = OUT_DIR / "summary.tsv"
VARIABLE_SITES_TSV = OUT_DIR / "variable_sites.tsv"
INDEX_MD = OUT_DIR / "index.md"


GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def read_fasta(path: Path) -> dict[str, str]:
    records = {}
    header = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def selected_m_records() -> dict[str, str]:
    cache = {}
    selected = {}
    for output_name, spec in SEQUENCES.items():
        source = spec["source"]
        if source not in cache:
            cache[source] = read_fasta(source)
        matches = [seq for header, seq in cache[source].items() if spec["match"](header)]
        if len(matches) != 1:
            raise RuntimeError(f"expected one match for {output_name}, found {len(matches)}")
        selected[output_name] = matches[0]
    return selected


def translate_m_glycoprotein_orf(sequence: str) -> tuple[str, int, int]:
    sequence = sequence.replace("-", "")
    stop_codons = {"TAA", "TAG", "TGA"}
    candidates = []
    for frame in range(3):
        active_start = None
        for index in range(frame, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]
            if active_start is None and codon == "ATG":
                active_start = index
            if codon in stop_codons:
                if active_start is not None:
                    orf_nt = index - active_start + 3
                    if orf_nt >= MIN_M_GLYCOPROTEIN_ORF_NT:
                        candidates.append((orf_nt, active_start, index + 3))
                active_start = None
    if not candidates:
        raise RuntimeError("could not infer a long M glycoprotein ORF")

    _, start, end = max(candidates)
    cds = sequence[start:end]
    amino_acids = []
    for index in range(0, len(cds), 3):
        codon = cds[index : index + 3]
        amino_acids.append(GENETIC_CODE.get(codon, "X"))
    protein = "".join(amino_acids)
    if "*" in protein[:-1]:
        raise RuntimeError("internal stop codon found in translated M glycoprotein")
    return protein.rstrip("*"), start + 1, end


def run_mafft(input_fasta: Path, output_fasta: Path) -> None:
    with output_fasta.open("w", encoding="utf-8") as output:
        subprocess.run(["mafft", "--auto", str(input_fasta)], check=True, stdout=output)


def position_map(sequence: str) -> list[int | None]:
    pos = 0
    mapped = []
    for char in sequence:
        if char != "-":
            pos += 1
            mapped.append(pos)
        else:
            mapped.append(None)
    return mapped


def pairwise_summary(records: dict[str, str]) -> list[dict[str, str]]:
    rows = []
    for left, right in combinations(records, 2):
        compared = differences = gap_columns = ambiguous = 0
        for left_char, right_char in zip(records[left], records[right]):
            if left_char == "-" or right_char == "-":
                gap_columns += 1
                continue
            compared += 1
            if left_char not in "ACGT" and right_char not in "ACDEFGHIKLMNPQRSTVWY*":
                ambiguous += 1
            if right_char not in "ACGT" and right_char not in "ACDEFGHIKLMNPQRSTVWY*":
                ambiguous += 1
            if left_char != right_char:
                differences += 1
        identity = ((compared - differences) / compared * 100) if compared else 0
        rows.append(
            {
                "left": left,
                "right": right,
                "compared": str(compared),
                "differences": str(differences),
                "identity_percent": f"{identity:.3f}",
                "gap_columns": str(gap_columns),
                "ambiguous_compared": str(ambiguous),
            }
        )
    return rows


def write_summary(path: Path, nt_records: dict[str, str], aa_records: dict[str, str]) -> None:
    rows = []
    for kind, records in [("nucleotide", nt_records), ("protein", aa_records)]:
        for row in pairwise_summary(records):
            rows.append({"type": kind, **row})
    with path.open("w", encoding="utf-8") as handle:
        handle.write("type\tleft\tright\tcompared\tdifferences\tidentity_percent\tgap_columns\tambiguous_compared\n")
        for row in rows:
            handle.write(
                f"{row['type']}\t{row['left']}\t{row['right']}\t{row['compared']}\t"
                f"{row['differences']}\t{row['identity_percent']}\t{row['gap_columns']}\t"
                f"{row['ambiguous_compared']}\n"
            )


def write_variable_sites(path: Path, nt_records: dict[str, str], aa_records: dict[str, str]) -> None:
    names = list(SEQUENCES)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("type\talignment_column\tposition\tSwiss_2026_M\tNRC-4_2018_M_MN258192.1\tNRC-3_1997_M_MN258191.1\n")
        for kind, records in [("nucleotide", nt_records), ("protein", aa_records)]:
            maps = {name: position_map(records[name]) for name in names}
            for index, chars in enumerate(zip(*(records[name] for name in names)), start=1):
                if len(set(chars)) == 1:
                    continue
                position = next((maps[name][index - 1] for name in names if maps[name][index - 1]), "")
                handle.write(
                    f"{kind}\t{index}\t{position}\t"
                    f"{records[names[0]][index - 1]}\t{records[names[1]][index - 1]}\t{records[names[2]][index - 1]}\n"
                )


def load_fonts():
    try:
        return {
            "title": ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 24),
            "bold": ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18),
            "normal": ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 16),
            "small": ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12),
            "tiny": ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 10),
        }
    except OSError:
        fallback = ImageFont.load_default()
        return {name: fallback for name in ["title", "bold", "normal", "small", "tiny"]}


def draw_variable_site_plot(
    records: dict[str, str],
    output_png: Path,
    title: str,
    subtitle: str,
    position_label: str,
    chunk_size: int,
) -> None:
    names = list(SEQUENCES)
    labels = {name: SEQUENCES[name]["label"] for name in names}
    fonts = load_fonts()
    maps = {name: position_map(records[name]) for name in names}
    variable_columns = [
        index
        for index, chars in enumerate(zip(*(records[name] for name in names)))
        if len(set(chars)) > 1
    ]

    colors = {
        "A": (82, 142, 194),
        "C": (101, 174, 150),
        "G": (220, 171, 60),
        "T": (224, 105, 47),
        "-": (210, 210, 210),
        "default": (158, 188, 99),
    }
    text_color = (31, 38, 46)
    muted = (82, 92, 102)
    line = (205, 212, 218)
    white = (255, 255, 255)

    cell_w = 18 if chunk_size > 40 else 28
    cell_h = 26
    label_w = 140
    left = 34
    top = 34
    chunks = [variable_columns[i : i + chunk_size] for i in range(0, len(variable_columns), chunk_size)] or [[]]
    width = max(
        left + label_w + max(len(chunk) for chunk in chunks) * cell_w + 44,
        920,
    )
    height = top + 140 + len(chunks) * (len(names) * cell_h + 86) + 76

    image = Image.new("RGB", (width, height), white)
    draw = ImageDraw.Draw(image)
    draw.text((left, top), title, fill=text_color, font=fonts["title"])
    draw.text((left, top + 34), subtitle, fill=muted, font=fonts["small"])

    ruler_x = left + label_w
    ruler_y = top + 88
    ruler_w = width - ruler_x - 40
    max_pos = max(pos for mapped in maps.values() for pos in mapped if pos)
    draw.text((left, ruler_y - 8), position_label, fill=text_color, font=fonts["small"])
    draw.line((ruler_x, ruler_y, ruler_x + ruler_w, ruler_y), fill=line, width=2)
    tick_step = 500 if max_pos > 1200 else 100
    for tick in range(0, max_pos + 1, tick_step):
        x = ruler_x + int((tick / max_pos) * ruler_w)
        draw.line((x, ruler_y - 6, x, ruler_y + 6), fill=(120, 128, 136), width=1)
        draw.text((x - 12, ruler_y + 10), str(tick), fill=muted, font=fonts["tiny"])
    for column in variable_columns:
        pos = next((maps[name][column] for name in names if maps[name][column]), None)
        if pos:
            x = ruler_x + int((pos / max_pos) * ruler_w)
            draw.line((x, ruler_y - 14, x, ruler_y - 8), fill=(199, 55, 45), width=1)

    legend_x = left
    legend_y = ruler_y + 34
    for char in ["A", "C", "G", "T", "-"]:
        draw.rectangle((legend_x, legend_y, legend_x + 16, legend_y + 16), fill=colors[char], outline=(120, 120, 120))
        draw.text((legend_x + 22, legend_y - 1), char, fill=text_color, font=fonts["small"])
        legend_x += 58
    draw.rectangle((legend_x, legend_y, legend_x + 16, legend_y + 16), fill=colors["default"], outline=(120, 120, 120))
    draw.text((legend_x + 22, legend_y - 1), "other amino acid", fill=text_color, font=fonts["small"])

    start_y = top + 158
    for chunk_index, columns in enumerate(chunks, start=1):
        if not columns:
            continue
        y = start_y + (chunk_index - 1) * (len(names) * cell_h + 86)
        draw.text(
            (left, y),
            f"Variable sites {columns[0] + 1}-{columns[-1] + 1} in alignment ({len(columns)} columns shown)",
            fill=text_color,
            font=fonts["bold"],
        )
        y += 30
        for cell_index, column in enumerate(columns):
            x = left + label_w + cell_index * cell_w
            pos = next((maps[name][column] for name in names if maps[name][column]), None)
            if cell_index % 4 == 0 or cell_index == len(columns) - 1:
                draw.text((x - 2, y), str(pos), fill=muted, font=fonts["tiny"])
        y += 18
        for name in names:
            draw.text((left, y + 4), labels[name], fill=text_color, font=fonts["normal"])
            for cell_index, column in enumerate(columns):
                x = left + label_w + cell_index * cell_w
                char = records[name][column]
                fill = colors.get(char, colors["default"])
                draw.rectangle((x, y, x + cell_w - 2, y + cell_h - 2), fill=fill, outline=white)
                draw.text((x + 4, y + 4), char, fill=text_color, font=fonts["small"])
            y += cell_h

    image.save(output_png)


def write_index() -> None:
    INDEX_MD.write_text(
        "\n".join(
            [
                "source: https://nextstrain.org/groups/hodcroftlab/andv/M",
                "method: MAFFT v7.505",
                "segment: M",
                "reason: M encodes the Andes virus glycoprotein precursor, relevant to entry/tropism comparisons.",
                "",
                "input:",
                "- Swiss_2026_M from data/fasta/ANDV-Switzerland-Hu-3337-2026/ANDV-Switzerland-Hu-3337-2026.fasta",
                "- NRC-4_2018_M_MN258192.1 from data/fasta/nextstrain-tips/M.fasta",
                "- NRC-3_1997_M_MN258191.1 from data/fasta/nextstrain-tips/M.fasta",
                "",
                "output:",
                "- input.fasta",
                "- aligned.fasta",
                "- M_closest_variable_sites.png",
                "- protein_input.fasta",
                "- protein_aligned.fasta",
                "- M_glycoprotein_closest_variable_sites.png",
                "- summary.tsv",
                "- variable_sites.tsv",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    nt_records = selected_m_records()
    write_fasta(NT_INPUT_FASTA, nt_records)
    run_mafft(NT_INPUT_FASTA, NT_ALIGNED_FASTA)
    nt_aligned = read_fasta(NT_ALIGNED_FASTA)

    inferred_orfs = {}
    aa_records = {}
    for name, sequence in nt_records.items():
        protein, start, end = translate_m_glycoprotein_orf(sequence)
        aa_records[name] = protein
        inferred_orfs[name] = (start, end, len(protein))
    write_fasta(AA_INPUT_FASTA, aa_records)
    run_mafft(AA_INPUT_FASTA, AA_ALIGNED_FASTA)
    aa_aligned = read_fasta(AA_ALIGNED_FASTA)

    write_summary(SUMMARY_TSV, nt_aligned, aa_aligned)
    write_variable_sites(VARIABLE_SITES_TSV, nt_aligned, aa_aligned)
    draw_variable_site_plot(
        nt_aligned,
        NT_VARIABLE_PNG,
        "ANDV M segment: Swiss 2026 vs closest NRC relatives",
        "Variable nucleotide alignment columns only. Grey marks gaps.",
        "M nt position",
        chunk_size=58,
    )
    draw_variable_site_plot(
        aa_aligned,
        AA_VARIABLE_PNG,
        "ANDV M glycoprotein: Swiss 2026 vs closest NRC relatives",
        "Variable amino-acid alignment columns only. Translation uses the inferred long M glycoprotein ORF.",
        "GPC aa position",
        chunk_size=40,
    )
    write_index()
    with INDEX_MD.open("a", encoding="utf-8") as handle:
        handle.write("inferred_m_glycoprotein_orfs:\n")
        for name, (start, end, protein_length) in inferred_orfs.items():
            handle.write(f"- {name}: nt {start}-{end}, {protein_length} aa\n")

    print(f"wrote {NT_ALIGNED_FASTA}")
    print(f"wrote {NT_VARIABLE_PNG}")
    print(f"wrote {AA_ALIGNED_FASTA}")
    print(f"wrote {AA_VARIABLE_PNG}")
    print(f"wrote {SUMMARY_TSV}")
    print(f"wrote {VARIABLE_SITES_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
