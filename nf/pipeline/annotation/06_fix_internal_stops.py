#!/usr/bin/env python3
"""Detect and repair CDS models with an in-frame stop codon in a funannotate/tbl2asn .tbl file.

Background: funannotate can emit a gene model whose spliced CDS has an internal stop codon,
often a stop formed across a splice junction (e.g. exon 1 ends in T, exon 2 starts with GA).
tbl2asn reports it as SEQ_FEAT.InternalStop and SEQ_INST.StopInProtein. The model is
internally consistent (canonical GT..AG, CDS length divisible by 3), so only the translation
shows the problem.

What this script does, per CDS with an internal stop:
  1. Try to move ONE splice site (donor or acceptor of one intron) by up to --window nt so
     that the CDS has an ATG start (if it had one), a single terminal stop, length divisible
     by 3, and a GT..AG intron at the moved site.
  2. Accept a move only when the new intron is supported by evidence introns (--introns,
     GFF3 exon features such as the PASA or transcript alignments) and the original intron
     is not. Exactly one such move must exist. Without evidence the sequence alone cannot
     say which move is right, so no guess is made.
  3. Otherwise drop the whole gene (gene, mRNA, CDS) from the .tbl, unless --no-drop is set.

Modes:
  --check    only report; exit 1 if any internal stop is found.
  (default)  write the repaired table to --out-tbl and a TSV report to --report.
Exit codes: 0 all resolved or none found, 1 unresolved genes remain, 2 tbl2asn (--val)
flagged a gene that this script did not detect (the tbl and the .val disagree).

Standard library only. Plain, .gz inputs are both accepted for --fasta and --introns.
Run with: /usr/bin/python3.12 06_fix_internal_stops.py ...

Typical use before tbl2asn (see 04_redo_sqn.sh):
  06_fix_internal_stops.py --tbl genome.tbl --fasta genome.fsa \
      --introns training/funannotate_train.pasa.gff3 training/funannotate_train.transcripts.gff3 \
      --out-tbl genome.fixed.tbl --report internal_stops.tsv
Then, after tbl2asn, gate on the validator output:
  06_fix_internal_stops.py --tbl genome.fixed.tbl --fasta genome.fsa --val genome.val --check
"""
import argparse
import collections
import gzip
import re
import sys

STOPS = {1: {"TAA", "TAG", "TGA"}, 4: {"TAA", "TAG"}}
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _c in enumerate(_AA):
    CODON[_B[_i // 16] + _B[(_i // 4) % 4] + _B[_i % 4]] = _c


def revcomp(s):
    return s.translate(COMP)[::-1]


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def num(x):
    return int(x.lstrip("<>"))


# ---------------------------------------------------------------- .tbl parsing

class Feature:
    def __init__(self, ftype, s, e):
        self.type = ftype
        self.ivs = [[s, e]]      # raw strings, may carry '<' or '>'
        self.quals = []          # raw qualifier lines

    def qual(self, key):
        pre = "\t\t\t" + key + "\t"
        return [q[len(pre):] for q in self.quals if q.startswith(pre)]

    def ints(self):
        return [(num(s), num(e)) for s, e in self.ivs]


def parse_tbl(path):
    """Return [{'name', 'groups'}]; a group is a gene and its features (the first group
    holds anything before the first gene, e.g. the REFERENCE line)."""
    contigs = []
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Feature"):
                cur = {"name": line.split()[1], "header": line, "groups": [[]]}
                contigs.append(cur)
            elif not line.strip():
                continue
            elif line.startswith("\t\t\t"):
                cur["groups"][-1][-1].quals.append(line)
            else:
                cols = line.split("\t")
                if len(cols) >= 3 and cols[2]:
                    feat = Feature(cols[2], cols[0], cols[1])
                    if feat.type == "gene":
                        cur["groups"].append([])
                    cur["groups"][-1].append(feat)
                else:
                    cur["groups"][-1][-1].ivs.append([cols[0], cols[1]])
    return contigs


def write_tbl(contigs):
    out = []
    for c in contigs:
        out.append(c["header"])
        for g in c["groups"]:
            for f in g:
                out.append("%s\t%s\t%s" % (f.ivs[0][0], f.ivs[0][1], f.type))
                out.extend("%s\t%s" % (s, e) for s, e in f.ivs[1:])
                out.extend(f.quals)
    return "\n".join(out) + "\n"


def read_fasta(path):
    seqs, cur = {}, None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


# ---------------------------------------------------------------- sequence logic

def strand_of(gene):
    s, e = gene.ints()[0]
    return "-" if s > e else "+"


def cds_seq(ivs, genome, strand):
    """Spliced CDS sequence for tx-ordered (start, end) intervals."""
    parts = []
    for s, e in ivs:
        parts.append(genome[s - 1:e] if s <= e else revcomp(genome[e - 1:s]))
    return "".join(parts)


def frame_of(feat):
    cs = feat.qual("codon_start")
    return int(cs[0]) if cs else 1


def internal_stops(seq, cs, stops):
    """List of (aa_position, codon, nt_offset_in_seq) for stops before the last codon."""
    body = seq[cs - 1:]
    hits = []
    n = len(body) // 3
    for i in range(n):
        codon = body[3 * i:3 * i + 3]
        if codon in stops and not (i == n - 1):
            hits.append((i + 1, codon, 3 * i + cs - 1))
    return hits


def junction_of(ivs, nt_off):
    """1-based junction index whose boundary falls inside the 3 nt starting at nt_off, else 0."""
    total = 0
    for j, (s, e) in enumerate(ivs[:-1]):
        total += abs(e - s) + 1
        if nt_off < total < nt_off + 3:
            return j + 1
    return 0


def intron_of(ivs, j, strand):
    (s1, e1), (s2, e2) = ivs[j], ivs[j + 1]
    return (e1 + 1, s2 - 1) if strand == "+" else (s2 + 1, e1 - 1)


def splice_ok(genome, lo, hi, strand):
    seq = genome[lo - 1:hi]
    if strand == "-":
        seq = revcomp(seq)
    return seq[:2] == "GT" and seq[-2:] == "AG"


def load_introns(paths):
    """{(contig, lo, hi): {strands}} from exon features grouped by Parent."""
    exons = collections.defaultdict(list)
    for p in paths:
        with open_text(p) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "exon":
                    continue
                m = re.search(r"Parent=([^;]+)", f[8])
                if m:
                    exons[(p, f[0], f[6], m.group(1))].append((int(f[3]), int(f[4])))
    introns = collections.defaultdict(set)
    for (_, contig, strand, _), ex in exons.items():
        ex.sort()
        for (a, b), (c, d) in zip(ex, ex[1:]):
            if c - b > 1:
                introns[(contig, b + 1, c - 1)].add(strand)
    return introns


def supported(introns, contig, intron, strand):
    st = introns.get((contig,) + intron)
    return bool(st) and (strand in st or "." in st)


def find_shift(ivs, genome, contig, strand, stops, cs, introns, window, min_intron):
    """Return (candidates, reason). candidates: list of (new_ivs, junction, side, delta)."""
    sgn = 1 if strand == "+" else -1
    old = cds_seq(ivs, genome, strand)
    atg = old[cs - 1:cs + 2] == "ATG"
    found = []
    for j in range(len(ivs) - 1):
        if supported(introns, contig, intron_of(ivs, j, strand), strand):
            continue    # the original intron has evidence: this junction is not the culprit
        for side in (0, 1):
            for d in range(-window, window + 1):
                if d == 0:
                    continue
                new = [list(x) for x in ivs]
                if side == 0:
                    new[j][1] += d * sgn
                else:
                    new[j + 1][0] += d * sgn
                if any((e - s) * sgn < 0 for s, e in new):
                    continue
                lo, hi = intron_of(new, j, strand)
                if hi - lo + 1 < min_intron or not splice_ok(genome, lo, hi, strand):
                    continue
                if not supported(introns, contig, (lo, hi), strand):
                    continue
                seq = cds_seq(new, genome, strand)
                body = seq[cs - 1:]
                if len(body) % 3 or body[-3:] not in stops or (atg and body[:3] != "ATG"):
                    continue
                if internal_stops(seq, cs, stops):
                    continue
                found.append(([tuple(x) for x in new], j + 1, "donor" if side == 0 else "acceptor", d))
    return found


def apply_shift(cds, mrna, old_ivs, new_ivs):
    """Rewrite CDS (and the matching mRNA boundary) in place. Return False if mRNA cannot match."""
    changes = [(i, k, old_ivs[i][k], new_ivs[i][k]) for i in range(len(old_ivs)) for k in (0, 1)
               if old_ivs[i][k] != new_ivs[i][k]]
    if mrna is not None:
        mi = mrna.ints()
        for _, k, o, _ in changes:
            if not any(iv[k] == o for iv in mi):
                return False
        for _, k, o, n in changes:
            idx = [x for x, iv in enumerate(mi) if iv[k] == o][0]
            mrna.ivs[idx][k] = str(n)
    for i, k, _, n in changes:
        cds.ivs[i][k] = str(n)
    return True


# ---------------------------------------------------------------- main

def parse_val(path, gene_ids):
    """Gene IDs that tbl2asn flagged with InternalStop/StopInProtein."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if "InternalStop" in line or "StopInProtein" in line:
                for m in re.finditer(r"gnl\|ncbi\|([^\]\s]+)\]", line):
                    ids.add(re.sub(r"-T\d+(_mrna)?$", "", m.group(1)))
                m = re.search(r"\((\S+) - ", line)
                if m:
                    ids.add(re.sub(r"-T\d+$", "", m.group(1)))
    return ids


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tbl", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-tbl")
    ap.add_argument("--report")
    ap.add_argument("--introns", nargs="*", default=[], help="GFF3 file(s) with exon evidence")
    ap.add_argument("--val", help="tbl2asn .val file to cross-check")
    ap.add_argument("--check", action="store_true", help="report only; exit 1 if any found")
    ap.add_argument("--no-drop", action="store_true", help="keep genes that cannot be repaired")
    ap.add_argument("--table", type=int, default=1, help="genetic code (1 or 4)")
    ap.add_argument("--window", type=int, default=30, help="max splice-site shift in nt")
    ap.add_argument("--min-intron", type=int, default=20)
    a = ap.parse_args(argv)
    if a.table not in STOPS:
        ap.error("--table must be one of %s" % sorted(STOPS))
    if not a.check and not a.out_tbl:
        ap.error("--out-tbl is required unless --check")
    stops = STOPS[a.table]

    contigs = parse_tbl(a.tbl)
    genome = read_fasta(a.fasta)
    introns = load_introns(a.introns) if a.introns else {}

    rows, unresolved, detected = [], 0, set()
    for c in contigs:
        keep = []
        for g in c["groups"]:
            gene = next((f for f in g if f.type == "gene"), None)
            if gene is None:
                keep.append(g)
                continue
            locus = (gene.qual("locus_tag") or ["?"])[0]
            mrnas = {(f.qual("transcript_id") or [None])[0]: f for f in g if f.type == "mRNA"}
            drop = False
            for cds in [f for f in g if f.type == "CDS"]:
                if c["name"] not in genome:
                    sys.exit("ERROR: contig %s is in the tbl but not in %s" % (c["name"], a.fasta))
                strand = strand_of(gene)
                ivs = cds.ints()
                cs = frame_of(cds)
                seq = cds_seq(ivs, genome[c["name"]], strand)
                hits = internal_stops(seq, cs, stops)
                if not hits:
                    continue
                detected.add(locus)
                tid = (cds.qual("transcript_id") or ["?"])[0]
                aa, codon, off = hits[0]
                row = [locus, tid, c["name"], strand, str(aa), codon,
                       str(junction_of(ivs, off) or "within_exon"), "", ""]
                if a.check:
                    row[7] = "detected"
                    rows.append(row)
                    continue
                cands = find_shift(ivs, genome[c["name"]], c["name"], strand, stops, cs,
                                   introns, a.window, a.min_intron) if len(ivs) > 1 else []
                if len(hits) > 1:
                    cands = []
                    row[8] = "%d internal stops; not repaired" % len(hits)
                if len(cands) == 1 and apply_shift(cds, mrnas.get(tid), ivs, [list(x) for x in cands[0][0]]):
                    new_ivs, j, side, d = cands[0]
                    row[7], row[8] = "fixed", "junction %d %s moved %+d nt (intron supported by evidence)" % (j, side, d)
                elif a.no_drop:
                    row[7] = "unresolved"
                    row[8] = row[8] or ("%d candidate(s) with evidence" % len(cands))
                    unresolved += 1
                else:
                    row[7] = "dropped"
                    row[8] = row[8] or ("%d candidate(s) with evidence; gene removed" % len(cands))
                    drop = True
                rows.append(row)
            if not drop:
                keep.append(g)
        c["groups"] = keep

    hdr = ["gene", "transcript", "contig", "strand", "first_stop_aa", "codon", "junction", "action", "detail"]
    if a.report:
        with open(a.report, "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            fh.writelines("\t".join(r) + "\n" for r in rows)
    for r in rows:
        print("\t".join(r))
    n_fixed = sum(r[7] == "fixed" for r in rows)
    n_drop = sum(r[7] == "dropped" for r in rows)
    print("internal-stop genes: %d (fixed %d, dropped %d, unresolved %d)" %
          (len(rows), n_fixed, n_drop, unresolved + (len(rows) if a.check else 0)), file=sys.stderr)

    if a.out_tbl and not a.check:
        with open(a.out_tbl, "w") as fh:
            fh.write(write_tbl(contigs))

    if a.val:
        missed = parse_val(a.val, detected) - detected
        if missed:
            print("ERROR: tbl2asn flagged genes this script did not detect: %s" %
                  ", ".join(sorted(missed)), file=sys.stderr)
            return 2
    if a.check:
        return 1 if rows else 0
    return 1 if unresolved else 0


if __name__ == "__main__":
    sys.exit(main())
