# =========================
# BioGen Mini-Project — biogen_gui.py
# Modern colorful Tkinter GUI calling main.py
# =========================

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import main  # algorithms file


# --------------------------
# Helpers
# --------------------------
def safe_set_text(widget: tk.Text, text: str) -> None:
    widget.config(state="normal")
    widget.delete("1.0", "end")
    widget.insert("1.0", text)
    widget.config(state="disabled")


def get_text(widget: tk.Text) -> str:
    return widget.get("1.0", "end").strip()


def extract_seq_from_text(widget: tk.Text) -> str:
    raw = get_text(widget)
    seq_lines = []
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        seq_lines.append(line)
    return "".join(seq_lines)


# --------------------------
# GUI App
# --------------------------
class BioGenGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("BioGen — DNA Analysis & Mutation Detection")
        self.geometry("1150x740")
        self.minsize(1050, 680)

        # Data storage for batch reports
        self.batch_reports = []

        # Apply colorful modern style
        self.apply_theme()

        # Notebook (tabs)
        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True, padx=12, pady=12)

        self.tab_single = ttk.Frame(notebook, padding=12)
        self.tab_compare = ttk.Frame(notebook, padding=12)
        self.tab_batch = ttk.Frame(notebook, padding=12)
        self.tab_3d = ttk.Frame(notebook, padding=12)

        notebook.add(self.tab_single, text=" Single ")
        notebook.add(self.tab_compare, text=" Compare ")
        notebook.add(self.tab_batch, text=" FASTA Batch ")
        notebook.add(self.tab_3d, text=" 3D Viewer ")

        self.build_single_tab()
        self.build_compare_tab()
        self.build_batch_tab()
        self.build_3d_tab()

    # --------------------------
    # Modern colorful theme (ttk)
    # --------------------------
    def apply_theme(self):
        self.configure(bg="#0b1220")  # deep navy

        style = ttk.Style(self)
        style.theme_use("clam")

        # Palette
        BG = "#0b1220"
        CARD = "#101a32"
        CARD2 = "#0f1a2e"
        TEXT = "#e8ecff"
        MUTED = "#aab3d6"
        ACCENT = "#7c5cff"  # purple
        ACCENT2 = "#00d4ff" # cyan
        OK = "#2dd4bf"      # teal
        WARN = "#ffb020"    # orange

        # Base styles
        style.configure(".", background=BG, foreground=TEXT, fieldbackground=CARD, bordercolor=CARD)
        style.configure("TFrame", background=BG)
        style.configure("Card.TFrame", background=CARD)
        style.configure("Card2.TFrame", background=CARD2)

        style.configure("TLabel", background=BG, foreground=TEXT, font=("Segoe UI", 10))
        style.configure("Title.TLabel", font=("Segoe UI", 14, "bold"), foreground=TEXT)
        style.configure("Sub.TLabel", font=("Segoe UI", 10), foreground=MUTED)
        style.configure("Badge.TLabel", font=("Segoe UI", 9, "bold"), foreground=BG, background=OK)

        style.configure("TButton", font=("Segoe UI", 10, "bold"), padding=(12, 8), background=CARD, foreground=TEXT)
        style.map(
            "TButton",
            background=[("active", "#18264a")],
            foreground=[("active", TEXT)]
        )

        style.configure("Accent.TButton", background=ACCENT, foreground="white")
        style.map("Accent.TButton", background=[("active", "#6a4cff")])

        style.configure("Cyan.TButton", background=ACCENT2, foreground=BG)
        style.map("Cyan.TButton", background=[("active", "#00c2ea")])

        # Notebook tabs
        style.configure("TNotebook", background=BG, borderwidth=0)
        style.configure("TNotebook.Tab", background=CARD, foreground=TEXT, padding=(14, 10), font=("Segoe UI", 10, "bold"))
        style.map("TNotebook.Tab", background=[("selected", ACCENT)], foreground=[("selected", "white")])

        # Treeview (batch list)
        style.configure("Treeview",
                        background=CARD, fieldbackground=CARD,
                        foreground=TEXT, rowheight=28, bordercolor=CARD, borderwidth=0,
                        font=("Segoe UI", 10))
        style.configure("Treeview.Heading",
                        background="#18264a",
                        foreground=TEXT,
                        font=("Segoe UI", 10, "bold"))
        style.map("Treeview", background=[("selected", "#263a7a")])

        # Save palette for later usage if needed
        self._ui = {"BG": BG, "CARD": CARD, "TEXT": TEXT, "MUTED": MUTED, "ACCENT": ACCENT, "ACCENT2": ACCENT2, "OK": OK, "WARN": WARN}

    # --------------------------
    # TAB 1 — Single
    # --------------------------
    def build_single_tab(self):
        header = ttk.Frame(self.tab_single, style="TFrame")
        header.pack(fill="x", pady=(0, 10))

        ttk.Label(header, text="Single DNA Analysis", style="Title.TLabel").pack(anchor="w")
        ttk.Label(header, text="Paste a DNA sequence or load a FASTA (first record).", style="Sub.TLabel").pack(anchor="w", pady=(2, 0))

        body = ttk.Frame(self.tab_single, style="TFrame")
        body.pack(fill="both", expand=True)

        left = ttk.Frame(body, style="Card.TFrame", padding=12)
        right = ttk.Frame(body, style="Card.TFrame", padding=12)
        left.pack(side="left", fill="both", expand=True, padx=(0, 10))
        right.pack(side="right", fill="both", expand=True)

        ttk.Label(left, text="Input", style="Title.TLabel").pack(anchor="w")
        self.single_input = tk.Text(left, height=18, wrap="word", bg=self._ui["CARD"], fg=self._ui["TEXT"], insertbackground="white", relief="flat")
        self.single_input.pack(fill="both", expand=True, pady=(8, 10))

        btns = ttk.Frame(left, style="Card.TFrame")
        btns.pack(fill="x")
        ttk.Button(btns, text="Load FASTA (first record)", style="Cyan.TButton", command=self.load_single_fasta).pack(side="left")
        ttk.Button(btns, text="Analyze", style="Accent.TButton", command=self.run_single_analysis).pack(side="left", padx=8)
        ttk.Button(btns, text="Clear", command=lambda: self.single_input.delete("1.0", "end")).pack(side="left")

        ttk.Label(right, text="Results", style="Title.TLabel").pack(anchor="w")
        self.single_output = tk.Text(right, wrap="word", state="disabled",
                                     bg=self._ui["CARD"], fg=self._ui["TEXT"], relief="flat")
        self.single_output.pack(fill="both", expand=True, pady=(8, 0))

    def load_single_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose a FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            records = main.read_fasta_records(path)
            seq_id, seq = records[0]  # first record only
            self.single_input.delete("1.0", "end")
            self.single_input.insert("1.0", f">{seq_id}\n{seq}\n")
        except Exception as e:
            messagebox.showerror("FASTA Error", str(e))

    def run_single_analysis(self):
        seq = extract_seq_from_text(self.single_input)
        try:
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

            out = []
            out.append("✅ VALIDATION: OK (A/T/C/G only)\n")

            out.append("=== BASIC STATISTICS ===")
            out.append(f"Length: {int(stats['length'])}")
            out.append(f"A: {int(stats['A'])}  T: {int(stats['T'])}  C: {int(stats['C'])}  G: {int(stats['G'])}")
            out.append(f"GC%: {stats['GC%']:.2f}   AT%: {stats['AT%']:.2f}\n")

            out.append("=== REVERSE COMPLEMENT ===")
            out.append(revc + "\n")

            out.append("=== START/STOP POSITIONS (nt, 0-based) ===")
            out.append(f"Start (ATG): {starts if starts else 'None'}")
            out.append(f"Stop (TAA/TAG/TGA): {stops if stops else 'None'}\n")

            out.append("=== ORFs (frames 0/1/2) ===")
            out.append(f"ORF count: {len(orfs)}")
            for i, o in enumerate(orfs, start=1):
                out.append(f"- ORF {i}: frame={o.frame} start={o.start} end={o.end} len={len(o.dna)}")
            out.append("")

            out.append("=== PROTEINS FROM ORFs ===")
            if proteins:
                for i, p in enumerate(proteins, start=1):
                    out.append(f"- Protein {i} (len={len(p)}): {p}")
            else:
                out.append("No ORFs => no proteins.")
            out.append("")

            out.append("=== TRANSLATION (WHOLE SEQ) ===")
            out.append(f"Frame 0: {prot0}")
            out.append(f"Frame 1: {prot1}")
            out.append(f"Frame 2: {prot2}")

            safe_set_text(self.single_output, "\n".join(out))
        except Exception as e:
            messagebox.showerror("Analysis Error", str(e))

    # --------------------------
    # TAB 2 — Compare
    # --------------------------
    def build_compare_tab(self):
        header = ttk.Frame(self.tab_compare)
        header.pack(fill="x", pady=(0, 10))
        ttk.Label(header, text="Two Sequences Comparison", style="Title.TLabel").pack(anchor="w")
        ttk.Label(header, text="Load or paste two sequences (FASTA first record).", style="Sub.TLabel").pack(anchor="w", pady=(2, 0))

        body = ttk.Frame(self.tab_compare)
        body.pack(fill="both", expand=True)

        top = ttk.Frame(body)
        top.pack(fill="both", expand=True)

        left = ttk.Frame(top, style="Card.TFrame", padding=12)
        right = ttk.Frame(top, style="Card.TFrame", padding=12)
        left.pack(side="left", fill="both", expand=True, padx=(0, 10))
        right.pack(side="right", fill="both", expand=True)

        ttk.Label(left, text="Sequence 1 (Reference)", style="Title.TLabel").pack(anchor="w")
        self.cmp_input1 = tk.Text(left, height=12, wrap="word", bg=self._ui["CARD"], fg=self._ui["TEXT"], insertbackground="white", relief="flat")
        self.cmp_input1.pack(fill="both", expand=True, pady=(8, 10))
        ttk.Button(left, text="Load FASTA 1", style="Cyan.TButton", command=lambda: self.load_compare_fasta(self.cmp_input1)).pack(anchor="w")

        ttk.Label(right, text="Sequence 2 (Mutated)", style="Title.TLabel").pack(anchor="w")
        self.cmp_input2 = tk.Text(right, height=12, wrap="word", bg=self._ui["CARD"], fg=self._ui["TEXT"], insertbackground="white", relief="flat")
        self.cmp_input2.pack(fill="both", expand=True, pady=(8, 10))
        ttk.Button(right, text="Load FASTA 2", style="Cyan.TButton", command=lambda: self.load_compare_fasta(self.cmp_input2)).pack(anchor="w")

        actions = ttk.Frame(body)
        actions.pack(fill="x", pady=10)
        ttk.Button(actions, text="Compare", style="Accent.TButton", command=self.run_comparison).pack(side="left")
        ttk.Button(actions, text="Clear", command=self.clear_compare).pack(side="left", padx=8)

        result_card = ttk.Frame(body, style="Card.TFrame", padding=12)
        result_card.pack(fill="both", expand=True)

        ttk.Label(result_card, text="Comparison Results", style="Title.TLabel").pack(anchor="w")
        self.compare_output = tk.Text(result_card, wrap="word", state="disabled",
                                      bg=self._ui["CARD"], fg=self._ui["TEXT"], relief="flat")
        self.compare_output.pack(fill="both", expand=True, pady=(8, 0))

    def load_compare_fasta(self, target_text: tk.Text):
        path = filedialog.askopenfilename(
            title="Choose a FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            records = main.read_fasta_records(path)
            seq_id, seq = records[0]
            target_text.delete("1.0", "end")
            target_text.insert("1.0", f">{seq_id}\n{seq}\n")
        except Exception as e:
            messagebox.showerror("FASTA Error", str(e))

    def clear_compare(self):
        self.cmp_input1.delete("1.0", "end")
        self.cmp_input2.delete("1.0", "end")
        safe_set_text(self.compare_output, "")

    def run_comparison(self):
        seq1 = extract_seq_from_text(self.cmp_input1)
        seq2 = extract_seq_from_text(self.cmp_input2)

        try:
            lengths = main.compare_lengths(seq1, seq2)
            a1, a2, muts = main.detect_mutations(seq1, seq2)
            classified = main.classify_mutations(muts)
            summary = main.mutation_summary(seq1, muts)
            orf_cmp = main.compare_orfs(seq1, seq2)
            prot_cmp = main.compare_proteins(seq1, seq2)
            impacts = main.impact_analysis(seq1, seq2)
            view = main.alignment_view(seq1, seq2)

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

            out.append("=== MUTATIONS ===")
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

            safe_set_text(self.compare_output, "\n".join(out))

        except Exception as e:
            messagebox.showerror("Comparison Error", str(e))

    # --------------------------
    # TAB 3 — FASTA Batch (each sequence separated)
    # --------------------------
    def build_batch_tab(self):
        header = ttk.Frame(self.tab_batch)
        header.pack(fill="x", pady=(0, 10))

        ttk.Label(header, text="Batch FASTA Analysis", style="Title.TLabel").pack(anchor="w")
        ttk.Label(header, text="Load a FASTA containing many sequences. Click a row to see its report separately.",
                  style="Sub.TLabel").pack(anchor="w", pady=(2, 8))

        actions = ttk.Frame(self.tab_batch)
        actions.pack(fill="x", pady=(0, 10))

        ttk.Button(actions, text="Load multi-FASTA", style="Accent.TButton", command=self.load_batch_fasta).pack(side="left")
        ttk.Button(actions, text="Clear", command=self.clear_batch).pack(side="left", padx=8)

        summary_card = ttk.Frame(self.tab_batch, style="Card.TFrame", padding=12)
        summary_card.pack(fill="x", pady=(0, 10))

        ttk.Label(summary_card, text="File Summary", style="Title.TLabel").pack(anchor="w")
        self.batch_summary = tk.Text(summary_card, height=4, wrap="word", state="disabled",
                                     bg=self._ui["CARD"], fg=self._ui["TEXT"], relief="flat")
        self.batch_summary.pack(fill="x", pady=(8, 0))

        mainbox = ttk.Frame(self.tab_batch)
        mainbox.pack(fill="both", expand=True)

        left = ttk.Frame(mainbox, style="Card.TFrame", padding=12)
        right = ttk.Frame(mainbox, style="Card.TFrame", padding=12)
        left.pack(side="left", fill="both", expand=False, padx=(0, 10))
        right.pack(side="right", fill="both", expand=True)

        ttk.Label(left, text="Sequences List", style="Title.TLabel").pack(anchor="w")

        # Treeview table
        columns = ("id", "len", "gc", "orf")
        self.seq_table = ttk.Treeview(left, columns=columns, show="headings", height=18)
        self.seq_table.heading("id", text="ID")
        self.seq_table.heading("len", text="Length")
        self.seq_table.heading("gc", text="GC%")
        self.seq_table.heading("orf", text="ORFs")

        self.seq_table.column("id", width=260, anchor="w")
        self.seq_table.column("len", width=70, anchor="center")
        self.seq_table.column("gc", width=70, anchor="center")
        self.seq_table.column("orf", width=60, anchor="center")

        self.seq_table.pack(fill="both", expand=True, pady=(8, 0))
        self.seq_table.bind("<<TreeviewSelect>>", self.on_select_batch_row)

        ttk.Label(right, text="Selected sequence report", style="Title.TLabel").pack(anchor="w")
        self.batch_output = tk.Text(right, wrap="word", state="disabled",
                                    bg=self._ui["CARD"], fg=self._ui["TEXT"], relief="flat")
        self.batch_output.pack(fill="both", expand=True, pady=(8, 0))

    def clear_batch(self):
        self.batch_reports = []
        safe_set_text(self.batch_summary, "")
        safe_set_text(self.batch_output, "")
        for item in self.seq_table.get_children():
            self.seq_table.delete(item)

    def load_batch_fasta(self):
        path = filedialog.askopenfilename(
            title="Choose a multi-FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.txt"), ("All files", "*.*")]
        )
        if not path:
            return

        try:
            self.batch_reports = main.analyze_fasta_records(path)
            summary = main.summarize_fasta_reports(self.batch_reports)

            summary_text = []
            summary_text.append(f"Sequences: {summary['n_sequences']}")
            summary_text.append(f"Length: min={summary['min_len']}  max={summary['max_len']}  avg={summary['avg_len']:.2f}")
            summary_text.append(f"GC% average: {summary['avg_gc']:.2f}")
            safe_set_text(self.batch_summary, "\n".join(summary_text))

            # fill table
            for item in self.seq_table.get_children():
                self.seq_table.delete(item)

            for i, r in enumerate(self.batch_reports, start=1):
                seq_id = r["id"]
                length = int(r["stats"]["length"])
                gc = float(r["stats"]["GC%"])
                orf_count = int(r["orf_count"])
                self.seq_table.insert("", "end", iid=str(i-1), values=(seq_id, length, f"{gc:.2f}", orf_count))

            # auto-select first
            if self.batch_reports:
                self.seq_table.selection_set("0")
                self.on_select_batch_row()

        except Exception as e:
            messagebox.showerror("Batch FASTA Error", str(e))

    def on_select_batch_row(self, event=None):
        if not self.batch_reports:
            return

        sel = self.seq_table.selection()
        if not sel:
            return

        idx = int(sel[0])  # iid holds index
        r = self.batch_reports[idx]

        stats = r["stats"]
        proteins = r["proteins"]

        out = []
        out.append("=" * 70)
        out.append(f"ID: {r['id']}")
        out.append("=" * 70)
        out.append("=== BASIC STATISTICS ===")
        out.append(f"Length: {int(stats['length'])}")
        out.append(f"A: {int(stats['A'])}  T: {int(stats['T'])}  C: {int(stats['C'])}  G: {int(stats['G'])}")
        out.append(f"GC%: {stats['GC%']:.2f}   AT%: {stats['AT%']:.2f}")
        out.append("")
        out.append("=== ORFs / PROTEINS ===")
        out.append(f"ORF count: {int(r['orf_count'])}")
        out.append(f"Proteins from ORFs: {len(proteins)}")

        if proteins:
            out.append("")
            out.append("Proteins (showing up to 5):")
            for i, p in enumerate(proteins[:5], start=1):
                out.append(f"- Protein {i} (len={len(p)}): {p}")
            if len(proteins) > 5:
                out.append(f"... (+{len(proteins) - 5} more)")
        else:
            out.append("No ORFs => no proteins.")

        safe_set_text(self.batch_output, "\n".join(out))

    # --------------------------
    # TAB 4 — 3D Viewer
    # --------------------------
    def build_3d_tab(self):
        header = ttk.Frame(self.tab_3d)
        header.pack(fill="x", pady=(0, 10))
        ttk.Label(header, text="3D Protein Viewer (Mol*)", style="Title.TLabel").pack(anchor="w")
        ttk.Label(header, text="Open structures by PDB ID or AlphaFold UniProt ID.", style="Sub.TLabel").pack(anchor="w", pady=(2, 0))

        card = ttk.Frame(self.tab_3d, style="Card.TFrame", padding=12)
        card.pack(fill="x", pady=(10, 10))

        form = ttk.Frame(card, style="Card.TFrame")
        form.pack(fill="x")

        ttk.Label(form, text="PDB ID:", style="TLabel").grid(row=0, column=0, sticky="w", pady=6)
        self.pdb_entry = ttk.Entry(form, width=24)
        self.pdb_entry.grid(row=0, column=1, sticky="w", padx=8, pady=6)

        ttk.Button(form, text="Open PDB", style="Cyan.TButton", command=self.on_open_pdb).grid(row=0, column=2, sticky="w", padx=8, pady=6)

        ttk.Label(form, text="UniProt ID:", style="TLabel").grid(row=1, column=0, sticky="w", pady=6)
        self.uniprot_entry = ttk.Entry(form, width=24)
        self.uniprot_entry.grid(row=1, column=1, sticky="w", padx=8, pady=6)

        ttk.Button(form, text="Open AlphaFold", style="Cyan.TButton", command=self.on_open_uniprot).grid(row=1, column=2, sticky="w", padx=8, pady=6)

        ttk.Label(card, text="Examples: PDB = 1CRN | UniProt = P69905", style="Sub.TLabel").pack(anchor="w", pady=(8, 0))

        # optional helper
        helper = ttk.Frame(self.tab_3d, style="Card.TFrame", padding=12)
        helper.pack(fill="both", expand=True)

        ttk.Label(helper, text="Optional: Protein difference positions", style="Title.TLabel").pack(anchor="w")

        grid = ttk.Frame(helper, style="Card.TFrame")
        grid.pack(fill="both", expand=True, pady=(8, 0))

        ttk.Label(grid, text="Protein 1:", style="TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(grid, text="Protein 2:", style="TLabel").grid(row=0, column=1, sticky="w")

        self.prot1 = tk.Text(grid, height=6, wrap="word", bg=self._ui["CARD"], fg=self._ui["TEXT"], insertbackground="white", relief="flat")
        self.prot2 = tk.Text(grid, height=6, wrap="word", bg=self._ui["CARD"], fg=self._ui["TEXT"], insertbackground="white", relief="flat")
        self.prot1.grid(row=1, column=0, padx=(0, 10), sticky="nsew")
        self.prot2.grid(row=1, column=1, sticky="nsew")

        grid.columnconfigure(0, weight=1)
        grid.columnconfigure(1, weight=1)

        ttk.Button(helper, text="Show changed residue positions", style="Accent.TButton", command=self.on_changed_residues).pack(anchor="w", pady=(10, 0))
        self.changed_label = ttk.Label(helper, text="", style="Sub.TLabel")
        self.changed_label.pack(anchor="w", pady=(6, 0))

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
        p1 = get_text(self.prot1).replace(" ", "").replace("\n", "")
        p2 = get_text(self.prot2).replace(" ", "").replace("\n", "")
        if not p1 or not p2:
            messagebox.showinfo("Info", "Please paste Protein 1 and Protein 2 first.")
            return

        changed = main.changed_residue_positions(p1, p2)
        if not changed:
            self.changed_label.config(text="No differences (proteins are identical).")
        else:
            show = changed[:100]
            more = " ..." if len(changed) > 100 else ""
            self.changed_label.config(text=f"Changed residues (1-based): {show}{more}")


if __name__ == "__main__":
    app = BioGenGUI()
    app.mainloop()
