# =========================
# BioGen Mini-Project (Updated Full Code)
# Part 1 + Part 2 + Part 3 + 3D Viewer helpers
# + FASTA multi-record analysis (Batch)
# =========================

from __future__ import annotations
from dataclasses import dataclass
import webbrowser
import re

# =========================
# PART 1 — Single DNA Sequence Analysis
# =========================

DNA_ALPHABET = set("AGCT")

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}

CODON_TABLE = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",
    # Leucine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine (Start)
    "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",

    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",

    # Tyrosine
    "TAT": "Y", "TAC": "Y",
    # Histidine
    "CAT": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAT": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic Acid
    "GAT": "D", "GAC": "D",
    # Glutamic Acid
    "GAA": "E", "GAG": "E",

    # Cysteine
    "TGT": "C", "TGC": "C",
    # Tryptophan
    "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",

    # Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*",
}


def clean_seq(seq: str) -> str:
    """Remove spaces/newlines and uppercase."""
    return seq.replace(" ", "").replace("\n", "").upper()


def validate_dna(seq: str) -> None:
    """Raise ValueError if sequence contains characters not in A,T,C,G."""
    s = clean_seq(seq)
    bad = [base for base in s if base not in DNA_ALPHABET]
    if bad:
        bad_unique = sorted(set(bad))
        raise ValueError(f"Invalid DNA sequence: only A,T,C,G allowed. Found: {bad_unique}")


# =========================
# NEW — FASTA MULTI-RECORD SUPPORT (Batch)
# =========================

def read_fasta_records(path: str) -> list[tuple[str, str]]:
    """
    Read a FASTA file with MANY sequences.
    Returns list of (seq_id, dna_sequence).
    Supports multi-line sequences.

    Example FASTA:
    >seq1
    ATGC...
    >seq2
    TTAA...
    """
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # save previous record
                if current_id is not None:
                    seq = clean_seq("".join(seq_parts))
                    if seq:
                        validate_dna(seq)
                        records.append((current_id, seq))

                # start new record
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    # save last record
    if current_id is not None:
        seq = clean_seq("".join(seq_parts))
        if seq:
            validate_dna(seq)
            records.append((current_id, seq))

    if not records:
        raise ValueError("No FASTA records found (file empty or invalid).")

    return records


def analyze_fasta_records(path: str) -> list[dict]:
    """
    Analyze ALL sequences in FASTA:
    - stats
    - ORFs
    - proteins from ORFs
    Returns list of dict reports (one per sequence).
    """
    records = read_fasta_records(path)
    reports: list[dict] = []

    for seq_id, seq in records:
        stats = dna_statistics(seq)
        orfs = find_orfs(seq)
        proteins = extract_proteins_from_orfs(orfs)

        reports.append({
            "id": seq_id,
            "sequence": seq,       # keep/remove depending on teacher preference
            "stats": stats,
            "orf_count": len(orfs),
            "proteins": proteins,
        })

    return reports


def summarize_fasta_reports(reports: list[dict]) -> dict:
    """Global summary for the whole FASTA file."""
    n = len(reports)
    lengths = [r["stats"]["length"] for r in reports]
    gcs = [r["stats"]["GC%"] for r in reports]

    return {
        "n_sequences": n,
        "min_len": min(lengths),
        "max_len": max(lengths),
        "avg_len": sum(lengths) / n,
        "avg_gc": sum(gcs) / n,
    }


# =========================
# Part 1 Functions
# =========================

def sequence_length(seq: str) -> int:
    """Length after cleaning."""
    return len(clean_seq(seq))


def dna_statistics(seq: str) -> dict[str, float]:
    """Length, counts, GC%, AT%."""
    validate_dna(seq)
    s = clean_seq(seq)

    length = len(s)
    A_count = s.count("A")
    T_count = s.count("T")
    C_count = s.count("C")
    G_count = s.count("G")

    gc = ((C_count + G_count) / length) * 100 if length > 0 else 0
    at = ((A_count + T_count) / length) * 100 if length > 0 else 0

    return {
        "length": length,
        "A": A_count,
        "T": T_count,
        "C": C_count,
        "G": G_count,
        "GC%": gc,
        "AT%": at,
    }


def reverse_complement(seq: str) -> str:
    """Reverse-complement DNA."""
    validate_dna(seq)
    s = clean_seq(seq)

    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complement = "".join(comp[base] for base in s)
    return complement[::-1]


def split_codons(seq: str, frame: int = 0) -> list[str]:
    """Split into codons starting at frame 0/1/2."""
    validate_dna(seq)
    s = clean_seq(seq)

    codons = []
    for i in range(frame, len(s) - 2, 3):
        codons.append(s[i:i + 3])
    return codons


def start_positions(seq: str) -> list[int]:
    """Return nucleotide positions of 'ATG' in the sequence (0-based)."""
    validate_dna(seq)
    s = clean_seq(seq)

    pos = []
    for i in range(0, len(s) - 2):
        if s[i:i + 3] == START_CODON:
            pos.append(i)
    return pos


def stop_positions(seq: str) -> list[int]:
    """Return nucleotide positions of stop codons in the sequence (0-based)."""
    validate_dna(seq)
    s = clean_seq(seq)

    pos = []
    for i in range(0, len(s) - 2):
        if s[i:i + 3] in STOP_CODONS:
            pos.append(i)
    return pos


@dataclass
class ORF:
    frame: int
    start: int   # nucleotide index (0-based)
    end: int     # nucleotide index (exclusive)
    dna: str     # ORF DNA including stop codon


def find_orfs(seq: str) -> list[ORF]:
    """
    Find ORFs on forward strand in frames 0/1/2.
    ORF = ATG ... (stop codon TAA/TAG/TGA)
    Non-overlap strategy: once one ORF found, continue after its stop.
    """
    validate_dna(seq)
    s = clean_seq(seq)
    orfs: list[ORF] = []

    for frame in (0, 1, 2):
        codons = split_codons(s, frame)
        i = 0
        while i < len(codons):
            if codons[i] == START_CODON:
                j = i + 1
                while j < len(codons):
                    if codons[j] in STOP_CODONS:
                        start_nt = frame + i * 3
                        end_nt = frame + (j + 1) * 3  # include stop codon
                        dna_orf = s[start_nt:end_nt]
                        orfs.append(ORF(frame, start_nt, end_nt, dna_orf))
                        i = j + 1
                        break
                    j += 1
                else:
                    i += 1
            else:
                i += 1

    return orfs


def translate_orf(orf_dna: str) -> str:
    """Translate ORF DNA until stop codon (stop not included in protein)."""
    validate_dna(orf_dna)
    s = clean_seq(orf_dna)

    protein = []
    for i in range(0, len(s) - 2, 3):
        codon = s[i:i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)

    return "".join(protein)


def extract_proteins_from_orfs(orfs: list[ORF]) -> list[str]:
    """Translate every ORF DNA into protein."""
    return [translate_orf(o.dna) for o in orfs]


def translate_dna(seq: str, frame: int = 0, stop_at_stop: bool = False) -> str:
    """
    Translate the WHOLE sequence starting from a frame (0/1/2).
    If stop_at_stop=True, stop when encountering first stop codon.
    """
    validate_dna(seq)
    s = clean_seq(seq)

    protein = []
    for i in range(frame, len(s) - 2, 3):
        codon = s[i:i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*" and stop_at_stop:
            break
        protein.append(aa)

    return "".join(protein)


# =========================
# PART 2 — Two DNA Sequence Comparison
# =========================

def validate_two_sequences(seq1: str, seq2: str) -> tuple[str, str]:
    """Clean + validate both DNA sequences and return cleaned versions."""
    s1 = clean_seq(seq1)
    s2 = clean_seq(seq2)
    validate_dna(s1)
    validate_dna(s2)
    return s1, s2


def compare_lengths(seq1: str, seq2: str) -> dict:
    """Length comparison report."""
    s1, s2 = validate_two_sequences(seq1, seq2)
    len1, len2 = len(s1), len(s2)
    return {
        "len_seq1": len1,
        "len_seq2": len2,
        "same_length": (len1 == len2),
        "difference": abs(len1 - len2),
        "which_is_longer": "seq1" if len1 > len2 else "seq2" if len2 > len1 else "equal",
    }


@dataclass
class Mutation:
    type: str  # "substitution" | "insertion" | "deletion"
    pos: int   # position on reference seq1 (0-based). insertion is "after pos"
    ref: str   # base in seq1, or '-'
    alt: str   # base in seq2, or '-'


def global_align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -1) -> tuple[str, str]:
    """Global alignment (Needleman–Wunsch). Returns aligned strings with '-' gaps."""
    s1, s2 = validate_two_sequences(seq1, seq2)
    n, m = len(s1), len(s2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[None] * (m + 1) for _ in range(n + 1)]  # "D", "U", "L"

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        tb[i][0] = "U"
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
        tb[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            best = max(diag, up, left)
            score[i][j] = best

            if best == diag:
                tb[i][j] = "D"
            elif best == up:
                tb[i][j] = "U"
            else:
                tb[i][j] = "L"

    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            move = tb[i][j]
        elif i > 0:
            move = "U"
        else:
            move = "L"

        if move == "D":
            a1.append(s1[i - 1]); a2.append(s2[j - 1])
            i -= 1; j -= 1
        elif move == "U":
            a1.append(s1[i - 1]); a2.append("-")
            i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1])
            j -= 1

    return "".join(reversed(a1)), "".join(reversed(a2))


def detect_mutations(seq1: str, seq2: str) -> tuple[str, str, list[Mutation]]:
    """Align then detect substitutions/insertions/deletions."""
    a1, a2 = global_align(seq1, seq2)
    muts: list[Mutation] = []

    pos_ref = -1
    for k in range(len(a1)):
        c1, c2 = a1[k], a2[k]

        if c1 != "-":
            pos_ref += 1

        if c1 == c2:
            continue

        if c1 != "-" and c2 != "-":
            muts.append(Mutation("substitution", pos_ref, c1, c2))
        elif c1 != "-" and c2 == "-":
            muts.append(Mutation("deletion", pos_ref, c1, "-"))
        else:
            muts.append(Mutation("insertion", pos_ref, "-", c2))

    return a1, a2, muts


def classify_mutations(muts: list[Mutation]) -> list[dict]:
    """Add a class field: substitution -> transition/transversion."""
    out = []
    for m in muts:
        if m.type == "substitution":
            if (m.ref in PURINES and m.alt in PURINES) or (m.ref in PYRIMIDINES and m.alt in PYRIMIDINES):
                cls = "transition"
            else:
                cls = "transversion"
        else:
            cls = m.type

        out.append({"type": m.type, "pos": m.pos, "ref": m.ref, "alt": m.alt, "class": cls})
    return out


def mutation_summary(seq1: str, muts: list[Mutation]) -> dict:
    """Count mutations + compute mutation rate = total_mutations / length(reference)."""
    s1 = clean_seq(seq1)
    n = len(s1)

    counts = {"substitution": 0, "insertion": 0, "deletion": 0}
    for m in muts:
        counts[m.type] += 1

    total = sum(counts.values())
    rate = (total / n) if n > 0 else 0

    return {
        "length_ref": n,
        "total_mutations": total,
        "counts": counts,
        "mutation_rate": rate,
    }


def orf_report(seq: str) -> dict:
    """ORF detection + proteins."""
    s = clean_seq(seq)
    validate_dna(s)

    orfs = find_orfs(s)
    proteins = [translate_orf(o.dna) for o in orfs]

    return {"orf_count": len(orfs), "orfs": orfs, "proteins": proteins}


def compare_orfs(seq1: str, seq2: str) -> dict:
    """Compare ORFs before/after mutation."""
    r1 = orf_report(seq1)
    r2 = orf_report(seq2)

    orf_dna_1 = [o.dna for o in r1["orfs"]]
    orf_dna_2 = [o.dna for o in r2["orfs"]]

    common = set(orf_dna_1) & set(orf_dna_2)

    return {
        "orf_count_before": r1["orf_count"],
        "orf_count_after": r2["orf_count"],
        "common_orfs": len(common),
        "lost_orfs": len(set(orf_dna_1) - set(orf_dna_2)),
        "new_orfs": len(set(orf_dna_2) - set(orf_dna_1)),
    }


def compare_proteins(seq1: str, seq2: str) -> dict:
    """Compare proteins extracted from ORFs before/after."""
    r1 = orf_report(seq1)
    r2 = orf_report(seq2)

    p1 = r1["proteins"]
    p2 = r2["proteins"]

    common = set(p1) & set(p2)

    return {
        "proteins_before": len(p1),
        "proteins_after": len(p2),
        "common_proteins": len(common),
        "lost_proteins": len(set(p1) - set(p2)),
        "new_proteins": len(set(p2) - set(p1)),
    }


def mutation_in_orf(pos: int, orf: ORF) -> bool:
    """True if reference position is inside ORF region."""
    return orf.start <= pos < orf.end


def codon_index_in_orf(pos: int, orf: ORF) -> int:
    """Codon number (0-based) inside ORF for this nucleotide position."""
    return (pos - orf.start) // 3


def impact_analysis(seq1: str, seq2: str) -> list[dict]:
    """
    Protein impact analysis based on reference ORFs (seq1):
      - substitution inside ORF -> silent/missense/nonsense
      - indel inside ORF -> frameshift (simple assumption: 1-base indel)
    """
    s1, s2 = validate_two_sequences(seq1, seq2)
    _, _, muts = detect_mutations(s1, s2)

    orfs_ref = find_orfs(s1)
    impacts: list[dict] = []

    for m in muts:
        if m.type in ("insertion", "deletion"):
            inside = False
            for o in orfs_ref:
                if mutation_in_orf(m.pos, o):
                    inside = True
                    impacts.append({
                        "mutation": m,
                        "impact": "frameshift",
                        "orf_start": o.start,
                        "orf_end": o.end,
                    })
                    break
            if not inside:
                impacts.append({"mutation": m, "impact": "noncoding_indel"})
            continue

        found_orf = None
        for o in orfs_ref:
            if mutation_in_orf(m.pos, o):
                found_orf = o
                break

        if found_orf is None:
            impacts.append({"mutation": m, "impact": "noncoding_substitution"})
            continue

        idx = codon_index_in_orf(m.pos, found_orf)
        codon_start = found_orf.start + idx * 3
        codon_ref = s1[codon_start:codon_start + 3]

        codon_list = list(codon_ref)
        offset = m.pos - codon_start
        codon_list[offset] = m.alt
        codon_mut = "".join(codon_list)

        aa_ref = CODON_TABLE.get(codon_ref, "X")
        aa_mut = CODON_TABLE.get(codon_mut, "X")

        if aa_mut == "*":
            impact = "nonsense"
        elif aa_ref == aa_mut:
            impact = "silent"
        else:
            impact = "missense"

        impacts.append({
            "mutation": m,
            "impact": impact,
            "orf_start": found_orf.start,
            "orf_end": found_orf.end,
            "codon_ref": codon_ref,
            "codon_mut": codon_mut,
            "aa_ref": aa_ref,
            "aa_mut": aa_mut,
        })

    return impacts


# =========================
# PART 3 — Visualization (text)
# =========================

def alignment_symbols(aligned1: str, aligned2: str) -> str:
    """'|' for match, '*' for mutation."""
    return "".join("|" if c1 == c2 else "*" for c1, c2 in zip(aligned1, aligned2))


def highlight_alignment(aligned1: str, aligned2: str) -> tuple[str, str]:
    """Highlight mutations with brackets []."""
    h1, h2 = [], []
    for c1, c2 in zip(aligned1, aligned2):
        if c1 == c2:
            h1.append(c1); h2.append(c2)
        else:
            h1.append(f"[{c1}]"); h2.append(f"[{c2}]")
    return "".join(h1), "".join(h2)


def protein_compare_view(p1: str, p2: str) -> tuple[str, str, str]:
    """Protein comparison line-by-line (pad shorter with '-' then use '|' / '*')."""
    n = max(len(p1), len(p2))
    a1 = p1.ljust(n, "-")
    a2 = p2.ljust(n, "-")
    mid = "".join("|" if x == y else "*" for x, y in zip(a1, a2))
    return a1, mid, a2


def print_protein_view(p1: str, p2: str) -> None:
    """Console print for protein comparison."""
    a1, mid, a2 = protein_compare_view(p1, p2)
    print(a1)
    print(mid)
    print(a2)


def alignment_view(seq1: str, seq2: str) -> dict:
    """Returns alignment text block + highlighted."""
    a1, a2 = global_align(seq1, seq2)
    mid = alignment_symbols(a1, a2)
    h1, h2 = highlight_alignment(a1, a2)
    return {
        "aligned1": a1,
        "symbols": mid,
        "aligned2": a2,
        "highlighted1": h1,
        "highlighted2": h2,
    }


# =========================
# PART 3D — 3D Protein Visualization (Advanced Feature)
# =========================

def _is_valid_pdb_id(pdb_id: str) -> bool:
    """PDB IDs are often 4 alphanumeric chars (e.g., 1CRN)."""
    return bool(re.fullmatch(r"[0-9A-Za-z]{4}", pdb_id))


def _is_valid_uniprot_id(uniprot_id: str) -> bool:
    """Simple UniProt ID check (good enough for mini-project)."""
    return bool(re.fullmatch(r"[0-9A-Za-z]{6,10}", uniprot_id))


def open_3d_pdb(pdb_id: str) -> None:
    """Open Mol* viewer using a PDB ID."""
    pdb_id = pdb_id.strip().upper()
    if not pdb_id:
        raise ValueError("Please enter a PDB ID (example: 1CRN).")
    if not _is_valid_pdb_id(pdb_id):
        raise ValueError("Invalid PDB ID format. Example: 1CRN, 7VUI.")
    url = f"https://molstar.org/viewer/?pdb={pdb_id}"
    webbrowser.open(url)


def open_3d_alphafold(uniprot_id: str) -> None:
    """Open Mol* viewer using AlphaFold DB UniProt ID."""
    uniprot_id = uniprot_id.strip().upper()
    if not uniprot_id:
        raise ValueError("Please enter a UniProt ID (example: P69905).")
    if not _is_valid_uniprot_id(uniprot_id):
        raise ValueError("Invalid UniProt ID format. Example: P69905, Q8W3K0.")
    url = f"https://molstar.org/viewer/?afdb={uniprot_id}"
    webbrowser.open(url)


def changed_residue_positions(prot1: str, prot2: str) -> list[int]:
    """1-based residue indices where proteins differ (helps show mutation effects)."""
    n = min(len(prot1), len(prot2))
    changed = []

    for i in range(n):
        if prot1[i] != prot2[i]:
            changed.append(i + 1)

    for i in range(n, max(len(prot1), len(prot2))):
        changed.append(i + 1)

    return changed


# =========================
# Quick test (optional)
# =========================
if __name__ == "__main__":
    s1 = "ATGAAATTTTAA"
    s2 = "ATGAAAATTTTAA"  # insertion of A

    print("Part 1 stats:", dna_statistics(s1))
    print("Reverse complement:", reverse_complement(s1))
    print("ORFs:", find_orfs(s1))
    print("Translate frame0:", translate_dna(s1, 0))

    a1, a2, muts = detect_mutations(s1, s2)
    print("\nAligned:")
    print(a1)
    print(alignment_symbols(a1, a2))
    print(a2)
    print("\nMutations:", muts)
    print("Classified:", classify_mutations(muts))
    print("Summary:", mutation_summary(s1, muts))
    print("Impact:", impact_analysis(s1, s2))
