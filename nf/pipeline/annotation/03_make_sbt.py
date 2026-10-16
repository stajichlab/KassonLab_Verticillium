#!/usr/bin/env python3
"""Build per-BioSample SBT files from a template SBT plus a contact source SBT.

Taken from the TEMPLATE (e.g. Kasson_Vert.sbt): citation authors and their
affiliation, the paper title (Seqdesc pub), the BioProject and all other blocks.
Taken from the CONTACT source (an existing SAMN*.sbt): the Submit-block contact
(name and affiliation), the "ALT EMAIL" comment, and the BioSample entry (which
is replaced with each sample's own BioSample).

Writes <out-dir>/<BioSample>.sbt for every BioSample in the samples file.
The original files are copied once to <out-dir>/orig/ before they are replaced.

Usage (from the project root):
    pipeline/annotation/03_make_sbt.py --template lib/sbt/Kasson_Vert.sbt \
        --contact-from lib/sbt/orig/SAMN64094693.sbt --samples samples.csv --out-dir lib/sbt
"""
import argparse
import csv
import os
import re
import shutil
import sys

CIT = "  cit {\n"
DBLINK = re.compile(r'Seqdesc ::= user \{\n  type str "DBLink",.*?\n\}\n', re.S)
ALT = re.compile(r'data str "ALT EMAIL:[^"]*"')
BIOSAMPLE = re.compile(r'\n    \{\n      label str "BioSample",.*?\n    \}', re.S)


def split_contact(text, name):
    i = text.find(CIT)
    if not text.startswith("Submit-block ::= {\n") or i < 0:
        sys.exit(f"{name}: not the expected SBT layout (no Submit-block / cit block)")
    return text[:i], text[i:]


def build(template, contact_src, biosample):
    contact_part, _ = split_contact(contact_src, "contact source")
    _, tmpl_rest = split_contact(template, "template")
    out = contact_part + tmpl_rest
    m = DBLINK.search(contact_src)
    if not m or not DBLINK.search(out):
        sys.exit("DBLink block not found in contact source or template")
    block = BIOSAMPLE.sub(lambda _: "\n    {\n      label str \"BioSample\",\n      num 1,\n"
                          f"      data strs {{\n        \"{biosample}\"\n      }}\n    }}", m.group(0))
    if '"BioSample"' not in block:
        sys.exit("BioSample entry not found in contact source DBLink block")
    out = DBLINK.sub(lambda _: block, out, count=1)
    alt = ALT.search(contact_src)
    if alt:
        out = ALT.sub(lambda _: alt.group(0), out, count=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--contact-from", required=True)
    ap.add_argument("--samples", required=True, help="CSV with a BioSample column")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    template = open(a.template).read()
    contact_src = open(a.contact_from).read()
    with open(a.samples, newline="") as fh:
        samples = [r["BioSample"].strip() for r in csv.DictReader(fh) if r.get("BioSample", "").strip()]
    orig = os.path.join(a.out_dir, "orig")
    os.makedirs(orig, exist_ok=True)
    for bs in samples:
        dest = os.path.join(a.out_dir, f"{bs}.sbt")
        bak = os.path.join(orig, f"{bs}.sbt")
        if os.path.exists(dest) and not os.path.exists(bak):
            shutil.copy2(dest, bak)
        with open(dest, "w") as fh:
            fh.write(build(template, contact_src, bs))
        print(f"wrote {dest}")


if __name__ == "__main__":
    main()
