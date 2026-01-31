# biogen_gui.py
# GUI ONLY (Option B): calls functions from main.py

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import main  # <-- your algorithms file (Part 1/2/3/3D)


# --------------------------
# FASTA helper (GUI side)
# --------------------------
def read_fasta_file(path: str) -> tuple[str, str]:
    """
    Reads a FASTA file and returns (header, sequence).
    Supports multi-line sequence.
    """
    header = ""
    seq_parts = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].strip()
                # if there are multiple records, we just read the first one (simple)
                continue
            else:
                seq_parts.append(line)

    seq = "".join(seq_parts)
    if not seq:
        raise ValueError("FASTA file has no sequence content.")
    return header, seq


def safe_set_text(widget: tk.Text, text: str) -> None:
    widget.config(state="normal")
    widget.delete("1.0", "end")
    widget.insert("1.0", text)
    widget.config(state="disabled")


def get_text(widget: tk.Text) -> str:
    return widget.get("1.0", "end").strip()


# --------------------------
# GUI App
# --------------------------
class BioGenGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("BioGen — DNA Analysis & Mutation Detection")
        self.geometry("1100x700")

        # Style (simple, clean)
        style = ttk.Style(self)
        style.theme_use("clam")

        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True)

        self.tab_single = ttk.Frame(notebook, padding=10)
        self.tab_compare = ttk.Frame(notebook, padding=10)
        self.tab_3d = ttk.Frame(notebook, padding=10)

        notebook.add(self.tab_single, text="Single Analysis")
        notebook.add(self.tab_compare, text="Comparison")
        notebook.add(self.tab_3d, text="3D Viewer")

        self.build_single_tab()
        self.build_compare_tab()
        self.build_3d_tab()

    # --------------------------
    # TAB 1 — Single analysis
    # --------------------------
    def build_single_tab(self):
        left = ttk.Frame(self.tab_single)
        right = ttk.Frame(self.tab_single)

        left.pack(side="left", fill="both", expand=True, padx=(0, 10))
        right.pack(side="right", fill="both", expand=True)

        ttk.Label(left, text="DNA Input (FASTA or raw sequence)", font=("Arial", 12, "bold")).pack(anchor="w")

        self.single_input = tk.Text(left, height=18, wrap="word")
        self.single_input.pack(fill="both", expand=True, pady=8)

        btns = ttk.Frame(left)
        btns.pack(fill="x", pady=5)

        ttk.Button(btns, text="Load FASTA", command=self.load_single_fasta).pack(side="left")
        ttk.Button(btns, text="Analyze", command=self.run_single_analysis).pack(side="left", padx=8)
        ttk.Button(btns, text="Clear", command=lambda: self.single_input.delete("1.0", "end")).pack(side="left")

        ttk.Label(right, text="Results", font=("Arial", 12, "bold")).pack(anchor="w")
        self.single_output = tk.Text(right, height=28, wrap="word", state="disabled")
        self.single_output.pack(fill="both", expand=True, pady=8)

    def load_single_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose a FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            header, seq = read_fasta_file(path)
            self.single_input.delete("1.0", "end")
            self.single_input.insert("1.0", f">{header}\n{seq}\n")
        except Exception as e:
            messagebox.showerror("FASTA Error", str(e))

    def run_single_analysis(self):
        raw = get_text(self.single_input)

        # If FASTA text, remove header lines and keep sequence
        seq_lines = []
        for line in raw.splitlines():
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
        seq = "".join(seq_lines)

        try:
            # Part 1 calls
            main.validate_dna(seq)
            stats = main.dna_statistics(seq)
            revc = main.reverse_complement(seq)

            starts = main.start_positions(seq)
            stops = main.stop_positions(seq)

            orfs = main.find_orfs(seq)
            proteins = main.extract_proteins_from_orfs(orfs)

            prot0 = main.translate_dna(seq, frame=0, stop_at_stop=False)
            prot1 = main.translate_dna(seq, frame=1, stop_at_stop=False)
            prot2 = main.translate_dna(seq, frame=2, stop_at_stop=False)

            # Build output
            out = []
            out.append("=== VALIDATION ===\nOK: sequence contains only A/T/C/G.\n")

            out.append("=== BASIC STATISTICS ===")
            out.append(f"Length: {int(stats['length'])}")
            out.append(f"A: {int(stats['A'])}  T: {int(stats['T'])}  C: {int(stats['C'])}  G: {int(stats['G'])}")
            out.append(f"GC%: {stats['GC%']:.2f}   AT%: {stats['AT%']:.2f}\n")

            out.append("=== REVERSE COMPLEMENT ===")
            out.append(revc + "\n")

            out.append("=== CODONS (Start/Stop positions) ===")
            out.append(f"Start codon ATG positions (nt): {starts if starts else 'None'}")
            out.append(f"Stop codon positions (nt): {stops if stops else 'None'}\n")

            out.append("=== ORFs (forward strand, frames 0/1/2) ===")
            out.append(f"ORF count: {len(orfs)}")
            if orfs:
                for idx, o in enumerate(orfs, start=1):
                    out.append(f"- ORF {idx}: frame={o.frame}, start={o.start}, end={o.end}, dna_len={len(o.dna)}")
            out.append("")

            out.append("=== PROTEINS FROM ORFs ===")
            if proteins:
                for idx, p in enumerate(proteins, start=1):
                    out.append(f"- Protein {idx} (len={len(p)}): {p}")
            else:
                out.append("No ORFs => no proteins.\n")

            out.append("=== TRANSLATION (WHOLE SEQ) ===")
            out.append(f"Frame 0: {prot0}")
            out.append(f"Frame 1: {prot1}")
            out.append(f"Frame 2: {prot2}")

            safe_set_text(self.single_output, "\n".join(out))

        except Exception as e:
            messagebox.showerror("Analysis Error", str(e))

    # --------------------------
    # TAB 2 — Comparison
    # --------------------------
    def build_compare_tab(self):
        top = ttk.Frame(self.tab_compare)
        top.pack(fill="both", expand=True)

        left = ttk.Frame(top)
        right = ttk.Frame(top)
        left.pack(side="left", fill="both", expand=True, padx=(0, 10))
        right.pack(side="right", fill="both", expand=True)

        ttk.Label(left, text="Sequence 1 (Reference)", font=("Arial", 11, "bold")).pack(anchor="w")
        self.cmp_input1 = tk.Text(left, height=14, wrap="word")
        self.cmp_input1.pack(fill="both", expand=True, pady=6)

        btns1 = ttk.Frame(left)
        btns1.pack(fill="x", pady=3)
        ttk.Button(btns1, text="Load FASTA 1", command=lambda: self.load_compare_fasta(self.cmp_input1)).pack(side="left")

        ttk.Label(right, text="Sequence 2 (Mutated)", font=("Arial", 11, "bold")).pack(anchor="w")
        self.cmp_input2 = tk.Text(right, height=14, wrap="word")
        self.cmp_input2.pack(fill="both", expand=True, pady=6)

        btns2 = ttk.Frame(right)
        btns2.pack(fill="x", pady=3)
        ttk.Button(btns2, text="Load FASTA 2", command=lambda: self.load_compare_fasta(self.cmp_input2)).pack(side="left")

        actions = ttk.Frame(self.tab_compare)
        actions.pack(fill="x", pady=10)

        ttk.Button(actions, text="Compare", command=self.run_comparison).pack(side="left")
        ttk.Button(actions, text="Clear", command=self.clear_compare).pack(side="left", padx=8)

        ttk.Label(self.tab_compare, text="Comparison Results", font=("Arial", 12, "bold")).pack(anchor="w")
        self.compare_output = tk.Text(self.tab_compare, height=18, wrap="word", state="disabled")
        self.compare_output.pack(fill="both", expand=True, pady=6)

    def load_compare_fasta(self, target_text: tk.Text):
        path = filedialog.askopenfilename(
            title="Choose a FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            header, seq = read_fasta_file(path)
            target_text.delete("1.0", "end")
            target_text.insert("1.0", f">{header}\n{seq}\n")
        except Exception as e:
            messagebox.showerror("FASTA Error", str(e))

    def clear_compare(self):
        self.cmp_input1.delete("1.0", "end")
        self.cmp_input2.delete("1.0", "end")
        safe_set_text(self.compare_output, "")

    def run_comparison(self):
        seq1 = self.extract_seq_from_text(self.cmp_input1)
        seq2 = self.extract_seq_from_text(self.cmp_input2)

        try:
            # Part 2 / Part 3 calls
            lengths = main.compare_lengths(seq1, seq2)

            a1, a2, muts = main.detect_mutations(seq1, seq2)
            classified = main.classify_mutations(muts)
            summary = main.mutation_summary(seq1, muts)

            orf_cmp = main.compare_orfs(seq1, seq2)
            prot_cmp = main.compare_proteins(seq1, seq2)
            impacts = main.impact_analysis(seq1, seq2)

            view = main.alignment_view(seq1, seq2)  # aligned + symbols + highlighted

            # Build output text
            out = []
            out.append("=== LENGTHS ===")
            out.append(str(lengths))
            out.append("")

            out.append("=== ALIGNMENT (| match, * mutation) ===")
            out.append(view["aligned1"])
            out.append(view["symbols"])
            out.append(view["aligned2"])
            out.append("")
            out.append("Highlighted:")
            out.append(view["highlighted1"])
            out.append(view["highlighted2"])
            out.append("")

            out.append("=== MUTATIONS LIST ===")
            if not muts:
                out.append("No mutations detected.")
            else:
                for x in classified:
                    out.append(f"- {x['class']:12}  pos={x['pos']}  {x['ref']} -> {x['alt']}")
            out.append("")

            out.append("=== MUTATION SUMMARY ===")
            out.append(str(summary))
            out.append("")

            out.append("=== ORF COMPARISON ===")
            out.append(str(orf_cmp))
            out.append("")

            out.append("=== PROTEIN COMPARISON ===")
            out.append(str(prot_cmp))
            out.append("")

            out.append("=== IMPACT ANALYSIS ===")
            if not impacts:
                out.append("No impacts (no mutations).")
            else:
                for item in impacts:
                    m = item["mutation"]
                    out.append(f"- {item['impact']:18}  pos={m.pos}  {m.ref}->{m.alt}")
                    if "aa_ref" in item:
                        out.append(f"    codon {item['codon_ref']} -> {item['codon_mut']} | AA {item['aa_ref']} -> {item['aa_mut']}")
            out.append("")

            safe_set_text(self.compare_output, "\n".join(out))

        except Exception as e:
            messagebox.showerror("Comparison Error", str(e))

    def extract_seq_from_text(self, widget: tk.Text) -> str:
        raw = get_text(widget)
        seq_lines = []
        for line in raw.splitlines():
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
        return "".join(seq_lines)

    # --------------------------
    # TAB 3 — 3D Viewer
    # --------------------------
    def build_3d_tab(self):
        ttk.Label(self.tab_3d, text="3D Protein Visualization (Mol*)", font=("Arial", 12, "bold")).pack(anchor="w")

        box = ttk.Frame(self.tab_3d)
        box.pack(fill="x", pady=15)

        # PDB ID
        ttk.Label(box, text="PDB ID (RCSB):").grid(row=0, column=0, sticky="w")
        self.pdb_entry = ttk.Entry(box, width=20)
        self.pdb_entry.grid(row=0, column=1, sticky="w", padx=8)

        ttk.Button(box, text="Open PDB in 3D", command=self.on_open_pdb).grid(row=0, column=2, sticky="w", padx=8)

        # UniProt ID
        ttk.Label(box, text="UniProt ID (AlphaFold DB):").grid(row=1, column=0, sticky="w", pady=(10, 0))
        self.uniprot_entry = ttk.Entry(box, width=20)
        self.uniprot_entry.grid(row=1, column=1, sticky="w", padx=8, pady=(10, 0))

        ttk.Button(box, text="Open AlphaFold 3D", command=self.on_open_uniprot).grid(row=1, column=2, sticky="w", padx=8, pady=(10, 0))

        ttk.Label(self.tab_3d, text="Examples: PDB = 1CRN | UniProt = P69905").pack(anchor="w")

        # Optional: protein difference helper
        ttk.Separator(self.tab_3d).pack(fill="x", pady=15)
        ttk.Label(self.tab_3d, text="Optional: Protein difference positions", font=("Arial", 11, "bold")).pack(anchor="w")

        pbox = ttk.Frame(self.tab_3d)
        pbox.pack(fill="both", expand=True, pady=8)

        ttk.Label(pbox, text="Protein 1:").grid(row=0, column=0, sticky="w")
        ttk.Label(pbox, text="Protein 2:").grid(row=0, column=1, sticky="w")

        self.prot1 = tk.Text(pbox, height=6, width=50, wrap="word")
        self.prot2 = tk.Text(pbox, height=6, width=50, wrap="word")
        self.prot1.grid(row=1, column=0, padx=(0, 10), sticky="nsew")
        self.prot2.grid(row=1, column=1, sticky="nsew")

        pbox.columnconfigure(0, weight=1)
        pbox.columnconfigure(1, weight=1)

        ttk.Button(self.tab_3d, text="Show changed residue positions", command=self.on_changed_residues).pack(anchor="w", pady=5)

        self.changed_label = ttk.Label(self.tab_3d, text="")
        self.changed_label.pack(anchor="w")

    def on_open_pdb(self):
        try:
            main.open_3d_pdb(self.pdb_entry.get())
        except Exception as e:
            messagebox.showerror("3D Viewer Error", str(e))

    def on_open_uniprot(self):
        try:
            main.open_3d_alphafold(self.uniprot_entry.get())
        except Exception as e:
            messagebox.showerror("3D Viewer Error", str(e))

    def on_changed_residues(self):
        p1 = get_text(self.prot1).replace(" ", "").replace("\n", "").strip()
        p2 = get_text(self.prot2).replace(" ", "").replace("\n", "").strip()
        if not p1 or not p2:
            messagebox.showinfo("Info", "Please paste Protein 1 and Protein 2 first.")
            return

        changed = main.changed_residue_positions(p1, p2)
        if not changed:
            self.changed_label.config(text="No differences (proteins are identical).")
        else:
            # show only first 100 positions to keep UI readable
            show = changed[:100]
            more = " ..." if len(changed) > 100 else ""
            self.changed_label.config(text=f"Changed residues (1-based): {show}{more}")


if __name__ == "__main__":
    app = BioGenGUI()
    app.mainloop()
