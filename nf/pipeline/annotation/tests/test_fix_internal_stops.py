"""Tests for 06_fix_internal_stops.py.

Run from the repository root:
    /usr/bin/python3.12 -m unittest discover -s nf/pipeline/annotation/tests
"""
import importlib.util
import os
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "fix_stops", os.path.join(HERE, "..", "06_fix_internal_stops.py"))
fs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fs)

# True gene (plus strand): ATG AAA TTT CCC | intron | GGG AAA TAA
UP = "GATTACA" * 3
E1 = "ATGAAATTTCCC"
INTRON = "GTAAGT" + "C" * 31 + "TAG"   # 40 nt, GT..AG; the last 3 nt are an in-frame stop
E2 = "GGGAAATAA"
DOWN = "ACGT" * 10
GENOME = UP + E1 + INTRON + E2 + DOWN
E1S = len(UP) + 1
E1E = len(UP) + len(E1)
INE = E1E + len(INTRON)
E2S = INE + 1
E2E = INE + len(E2)


def make_tbl(acc_shift):
    """Model whose acceptor is acc_shift nt too far upstream."""
    s2 = E2S - acc_shift
    return "\n".join([
        ">Feature c1",
        "1\t%d\tREFERENCE" % len(GENOME),
        "\t\t\tNCBI\t12345",
        "%d\t%d\tgene" % (E1S, E2E),
        "\t\t\tlocus_tag\tG1",
        "%d\t%d\tmRNA" % (E1S, E1E),
        "%d\t%d" % (s2, E2E),
        "\t\t\tproduct\thypothetical protein",
        "\t\t\ttranscript_id\tgnl|ncbi|G1-T1_mrna",
        "\t\t\tprotein_id\tgnl|ncbi|G1-T1",
        "%d\t%d\tCDS" % (E1S, E1E),
        "%d\t%d" % (s2, E2E),
        "\t\t\tcodon_start\t1",
        "\t\t\tproduct\thypothetical protein",
        "\t\t\ttranscript_id\tgnl|ncbi|G1-T1_mrna",
        "\t\t\tprotein_id\tgnl|ncbi|G1-T1",
        "",
    ])


class Case(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.fsa = os.path.join(self.d, "g.fsa")
        with open(self.fsa, "w") as fh:
            fh.write(">c1\n%s\n" % GENOME)
        self.evid = os.path.join(self.d, "e.gff3")
        with open(self.evid, "w") as fh:
            fh.write("##gff-version 3\n")
            for a, b in ((E1S, E1E), (E2S, E2E)):
                fh.write("c1\tpasa\texon\t%d\t%d\t.\t+\t.\tID=x;Parent=t1\n" % (a, b))

    def write_tbl(self, shift):
        p = os.path.join(self.d, "g.tbl")
        with open(p, "w") as fh:
            fh.write(make_tbl(shift))
        return p

    def run_main(self, shift, extra=()):
        out = os.path.join(self.d, "out.tbl")
        rep = os.path.join(self.d, "rep.tsv")
        rc = fs.main(["--tbl", self.write_tbl(shift), "--fasta", self.fsa,
                      "--out-tbl", out, "--report", rep] + list(extra))
        return rc, out, rep

    def test_roundtrip(self):
        p = self.write_tbl(0)
        with open(p) as fh:
            self.assertEqual(fs.write_tbl(fs.parse_tbl(p)), fh.read())

    def test_detects_internal_stop(self):
        rc = fs.main(["--tbl", self.write_tbl(3), "--fasta", self.fsa, "--check"])
        self.assertEqual(rc, 1)

    def test_clean_gene_passes(self):
        rc = fs.main(["--tbl", self.write_tbl(0), "--fasta", self.fsa, "--check"])
        self.assertEqual(rc, 0)

    def test_fix_with_evidence_restores_true_model(self):
        rc, out, rep = self.run_main(3, ["--introns", self.evid])
        self.assertEqual(rc, 0)
        with open(out) as fh:
            self.assertEqual(fh.read(), make_tbl(0))
        with open(rep) as fh:
            self.assertIn("fixed", fh.read())

    def test_drop_without_evidence(self):
        rc, out, rep = self.run_main(3)
        self.assertEqual(rc, 0)
        with open(out) as fh:
            self.assertNotIn("G1", fh.read())
        with open(rep) as fh:
            self.assertIn("dropped", fh.read())

    def test_no_drop_leaves_gene_and_fails(self):
        rc, out, rep = self.run_main(3, ["--no-drop"])
        self.assertEqual(rc, 1)
        with open(out) as fh:
            self.assertIn("G1", fh.read())


if __name__ == "__main__":
    unittest.main()
