"""
Sup. Fig. 3: residues that distinguish annelid r-opsin1 from r-opsin3 orthologs,
taken from the alignment used for the trees (no re-alignment). Mollusc rhodopsins
polarise states (ancestral vs r-opsin1-/r-opsin3-specific); Clade B shown unclassified.

Usage:
  python3 Diagnostic_residues.py ALIGNMENT.fa LABELS.csv OUT_PREFIX [FULL_LENGTH.fa]

FULL_LENGTH.fa (optional, e.g. ../listA_filtered_add.fa) is only used to number
columns by full-length Platynereis r-opsin1 positions instead of 7TM-block positions.
"""
import csv, re, sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

matplotlib.rcParams["svg.fonttype"] = "none"      # editable text in Inkscape
matplotlib.rcParams["font.family"] = "DejaVu Sans"

aln_file, csv_file, out_prefix = sys.argv[1:4]
full_file = sys.argv[4] if len(sys.argv) > 4 else None

# ---------- read data ----------
lab = {}
with open(csv_file, newline="") as f:
    for r in csv.DictReader(f):
        lab[r["Seq_ID"]] = r["Tree_labels"]

def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split("/")[0].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line.upper()
    return seqs

seqs = read_fasta(aln_file)

def pick(keys):
    out = []
    for k in keys:                       # keep the order given here
        out += [s for s in seqs if k in lab.get(s, s) and s not in out]
    return out

R1  = pick(["Platynereis dumerilii r-opsin1", "Perionyx excavatus r-opsin1", "Malacoceros fuliginosus r-opsin1"])
OW  = pick(["CAH1779377"])
R3  = pick(["Platynereis dumerilii r-opsin3", "Platynereis dumerilii r-opsin2", "Malacoceros fuliginosus r-opsin3",
            "Capitella teleta r-opsin1", "Perionyx excavatus r-opsin3", "Helobdella"])
OUT = pick(["rhodops"])
CB  = pick(["Malacoceros fuliginosus r-opsin4", "Platynereis dumerilii r-opsin4", "Platynereis dumerilii r-opsin5",
            "Capitella teleta r-opsin2", "Capitella teleta r-opsin3",
            "CAH1795029", "CAH1794772", "CAH1795196"])
assert len(OW) == 1 and R1 and R3 and OUT and CB, "group lookup failed - check Tree_labels"
print("Clade B rows:", "; ".join(lab.get(s, s) for s in CB))
ref = R1[0]                               # Platynereis r-opsin1 = numbering reference

def states(group, i):
    return {seqs[s][i] for s in group if seqs[s][i] != "-"}

# ---------- diagnostic columns (same rules as diag_residues.py) ----------
diag = []                                 # (col, r1_states, r3_states)
for i in range(len(seqs[ref])):
    s1, s3 = states(R1, i), states(R3, i)
    if not s1 or not s3:
        continue
    if s1 & s3:
        c1 = [seqs[s][i] for s in R3]; c3 = [seqs[s][i] for s in R1]
        f1 = len(s1) == 1 and c1.count(next(iter(s1))) <= 1
        f3 = len(s3) == 1 and c3.count(next(iter(s3))) <= 1
        if not (f1 or f3):
            continue
        s1x, s3x = (s1, s3 - s1) if f1 else (s1 - s3, s3)
    else:
        s1x, s3x = s1, s3
    diag.append((i, s1x, s3x))

# ---------- reference positions (K296, NPxxY) from Platynereis r-opsin1 ----------
ung2col = [i for i, c in enumerate(seqs[ref]) if c != "-"]
ung = seqs[ref].replace("-", "")
ref_cols = {}
m = list(re.finditer(r"NP..Y", ung))
if m:
    n_pos = m[-1].start()
    ref_cols[ung2col[n_pos]] = "N302"
    ref_cols[ung2col[n_pos + 4]] = "Y306"
    if n_pos - 6 >= 0 and ung[n_pos - 6] == "K":
        ref_cols[ung2col[n_pos - 6]] = "K296"

# ---------- column numbering: full-length Pd r-opsin1 if available ----------
pos_lab = {c: k + 1 for k, c in enumerate(ung2col)}
numbering = "7TM block"
if full_file:
    full = read_fasta(full_file).get(ref, "").replace("-", "")
    j, mapped = 0, {}
    for k, aa in enumerate(ung):          # block residues as an ordered subsequence of full length
        while j < len(full) and full[j] != aa:
            j += 1
        if j == len(full):
            mapped = None; break
        mapped[ung2col[k]] = j + 1; j += 1
    if mapped:
        pos_lab, numbering = mapped, "full-length"
    else:
        print("WARNING: could not map block to full-length sequence; using 7TM-block numbering")

cols = sorted(set([d[0] for d in diag]) | set(ref_cols))
dinfo = {d[0]: d for d in diag}

# ---------- cell classes ----------
C_R1, C_R3, C_ANC, C_UNQ = "#5aae61", "#4393c3", "#bdbdbd", "#e08214"

def cell_class(s, i):
    a = seqs[s][i]
    if a == "-" or s in CB:              # Clade B shown for reference, not classified
        return None
    if i not in dinfo:                   # reference sites: flag only Owenia's unique residue (D at N302)
        return "unq" if (s == OW[0] and a not in states(R1 + R3 + OUT, i)) else None
    _, s1x, s3x = dinfo[i]
    anc = states(OUT, i)
    if a in anc:
        return "anc"
    if a in s1x:
        return "r1"
    if a in s3x:
        return "r3"
    if s == OW[0] and a not in states(R1 + R3, i):
        return "unq"
    return None

fill = {"r1": C_R1, "r3": C_R3, "anc": C_ANC, "unq": C_UNQ}

# ---------- table output ----------
with open(out_prefix + "_table.tsv", "w") as f:
    f.write("aln_col\tPdu_r1_pos\tref_site\tr1_states\tr3_states\tmollusc_states\tOwenia\tOwenia_class\n")
    for i in cols:
        d = dinfo.get(i)
        f.write(f"{i+1}\t{pos_lab.get(i,'')}\t{ref_cols.get(i,'')}\t"
                f"{''.join(sorted(d[1])) if d else ''}\t{''.join(sorted(d[2])) if d else ''}\t"
                f"{''.join(sorted(states(OUT,i)))}\t{seqs[OW[0]][i]}\t{cell_class(OW[0], i) or 'other'}\n")

# ---------- plot ----------
groups = [("r-opsin1", R1, "#E8F8DC"), ("", OW, "#E8F8DC"), ("r-opsin3", R3, "#D9FAFC"),
          ("Mollusca", OUT, "#F2F2F2"), ("Clade B", CB, "#F7F3FA")]
rows = []
for g, mem, bg in groups:
    rows += [(s, g, bg) for s in mem]
nR, nC = len(rows), len(cols)

cw, rh = 0.24, 0.24
fig_w = 3.6 + nC * cw
fig_h = 1.6 + nR * rh
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_xlim(-0.5, nC - 0.5)
ax.set_ylim(nR - 0.5, -2.2)
ax.axis("off")

for r, (s, g, bg) in enumerate(rows):
    ax.add_patch(Rectangle((-0.5, r - 0.5), nC, 1, facecolor=bg, edgecolor="none", zorder=0))
    for c, i in enumerate(cols):
        a = seqs[s][i]
        k = cell_class(s, i)
        if k:
            ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, facecolor=fill[k], alpha=0.65,
                                   edgecolor="white", linewidth=0.8, zorder=1))
        is_mfu = s.startswith("Mfu_")
        bold = s == OW[0] or is_mfu
        ax.text(c, r, a, ha="center", va="center", fontsize=7.5, family="DejaVu Sans Mono",
                fontweight="bold" if bold else "normal", color="#222222", zorder=2)
    name = lab.get(s, s).replace(",", "")
    ax.text(-0.8, r, name, ha="right", va="center", fontsize=7, style="italic",
            fontweight="bold" if bold else "normal",
            color="#b2182b" if is_mfu else "#222222")

# Owenia row outline
r_ow = [k for k, row in enumerate(rows) if row[0] == OW[0]][0]
ax.add_patch(Rectangle((-0.5, r_ow - 0.5), nC, 1, fill=False, edgecolor="#222222", linewidth=1.0, zorder=3))

# reference columns: dashed outline + label
for c, i in enumerate(cols):
    if i in ref_cols:
        ax.add_patch(Rectangle((c - 0.5, -0.5), 1, nR, fill=False, edgecolor="#555555",
                               linestyle=(0, (2, 1.5)), linewidth=0.9, zorder=3))
        ax.text(c, -1.75, ref_cols[i], ha="center", va="bottom", fontsize=6.5, rotation=90, color="#222222")
    ax.text(c, -0.65, str(pos_lab.get(i, "")), ha="center", va="bottom", fontsize=5.5, rotation=90, color="#555555")

# group bars on the right
start = 0
for g, mem, bg in groups:
    if mem and g:
        end = start + len(mem) - 1 + (len(OW) if g == "r-opsin1" else 0)
        ax.plot([nC - 0.2, nC - 0.2], [start - 0.35, end + 0.35], color="#222222", lw=1.6, clip_on=False)
        ax.text(nC + 0.05, (start + end) / 2, g, ha="left", va="center", fontsize=8, clip_on=False)
    start += len(mem)

nA = len(R1) + len(OW) + len(R3)
ax.plot([nC + 2.7, nC + 2.7], [-0.35, nA - 1 + 0.35], color="#222222", lw=1.6, clip_on=False)
ax.text(nC + 2.95, (nA - 1) / 2, "Clade A", ha="left", va="center", rotation=90, fontsize=9, clip_on=False)

handles = [Patch(facecolor=C_R1, alpha=0.65, label="r-opsin1-specific (derived)"),
           Patch(facecolor=C_R3, alpha=0.65, label="r-opsin3-specific (derived)"),
           Patch(facecolor=C_ANC, alpha=0.65, label="as in molluscs (ancestral)"),
           Patch(facecolor=C_UNQ, alpha=0.65, label="unique to O. fusiformis")]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.01), ncol=4, frameon=False,
          fontsize=7, handlelength=1.2, columnspacing=1.2)


for ext in ("png", "svg"):
    fig.savefig(f"{out_prefix}.{ext}", dpi=300, bbox_inches="tight")
print(f"{len(diag)} diagnostic columns + {len(ref_cols)} reference sites; numbering: {numbering}")
print(f"Owenia classes: " + ", ".join(f"{k}={sum(cell_class(OW[0], i) == k for i in dinfo)}"
                                      for k in ("r1", "r3", "anc", "unq")))
print("Saved:", out_prefix + ".png / .svg / _table.tsv")