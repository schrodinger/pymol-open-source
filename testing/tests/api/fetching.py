"""
fetch_by_uniprot, fetch_uniprot_ligands, list_uniprot_structures, fetch_fold
and color_plddt (modules/pymol/fetching.py). The servers are replaced by
canned answers; test_live asks the real ones and is not part of "--run all".
"""

import json
from unittest import mock

from pymol import CmdException, cmd, fetching, testing

XREF = {
    "uniProtKBCrossReferences": [
        {
            "database": "PDB",
            "id": "1TUP",
            "properties": [
                {"key": "Method", "value": "X-ray"},
                {"key": "Resolution", "value": "2.20 A"},
                {"key": "Chains", "value": "A/B/C=94-312"},
            ],
        },
        {"database": "PDB", "id": "2OCJ", "properties": []},
        {"database": "EMBL", "id": "X02469"},
    ]
}


def model(bfactors):
    """PDB text of a small peptide with the given B-factor per residue."""
    cmd.fab("A" * len(bfactors), "_model")
    for resi, b in enumerate(bfactors, 1):
        cmd.alter("_model and resi %d" % resi, "b=%r" % b)
    text = cmd.get_pdbstr("_model")
    cmd.delete("_model")
    return text


def colors(selection):
    found = set()
    cmd.iterate(selection, "found.add(color)", space={"found": found})
    return {cmd.get_color_tuple(c) for c in found}


class TestFetching(testing.PyMOLTestCase):
    def serve(self, answers):
        """Replace the HTTP call: `answers` maps a part of the URL to text."""

        def get(url, data=None):
            self.asked.append((url, data))
            return next((v for k, v in answers.items() if k in url), "")

        self.asked = []
        patch = mock.patch.object(fetching, "_get", get)
        patch.start()
        self.addCleanup(patch.stop)

    def test_uniprot(self):
        self.serve({"P04637.json": json.dumps(XREF)})
        self.assertEqual(cmd.list_uniprot_structures("p04637"), ["1TUP", "2OCJ"])
        self.assertEqual(
            fetching.pdb_entries("P04637")[0],
            ("1TUP", "X-ray", "2.20 A", "A/B/C=94-312"),
        )
        self.assertRaises(CmdException, cmd.list_uniprot_structures, "nonsense")

        def fetch(code, name, **_):
            if code == "2OCJ":
                raise CmdException("no such file")
            cmd.fragment("ala", name)

        with mock.patch.object(cmd, "fetch", fetch):
            self.assertEqual(cmd.fetch_by_uniprot("P04637"), ["P04637_1TUP"])
            self.assertEqual(cmd.get_object_list("P04637"), ["P04637_1TUP"])
            cmd.delete("*")
            names = cmd.fetch_by_uniprot("P04637", 1, "p53")
            self.assertEqual(cmd.get_object_list("p53"), names)
        self.serve({"P04637.json": json.dumps({})})
        self.assertRaises(CmdException, cmd.fetch_by_uniprot, "P04637")

    def test_uniprot_ligands(self):
        def component(name, weight):
            return {
                "nonpolymer_comp": {
                    "chem_comp": {"id": name, "name": name, "formula_weight": weight}
                }
            }

        cmd.fragment("benzene", "_lig")
        sdf = cmd.get_str("mol", "_lig")
        cmd.delete("_lig")
        entries = [
            {"nonpolymer_entities": [component("AZM", 222.2), component("ZN", 65.4)]},
            {"nonpolymer_entities": [component("BEZ", 122.1), component("XXX", 300)]},
            {"nonpolymer_entities": None},
        ]
        self.serve(
            {
                "P04637.json": json.dumps(XREF),
                "graphql": json.dumps({"data": {"entries": entries}}),
                "AZM_ideal": sdf,
                "BEZ_ideal": sdf,
                "ZN_ideal": sdf,
            }
        )
        # ZN is too light, XXX has no ideal model
        names = cmd.fetch_uniprot_ligands("p04637")
        self.assertEqual(names, ["P04637_AZM", "P04637_BEZ"])
        self.assertEqual(cmd.get_object_list("P04637_ligands"), names)
        self.assertEqual(cmd.count_atoms("P04637_AZM"), 12)
        # the two PDB ids of the accession were asked for in one request
        asked = [url for url, _ in self.asked if "graphql" in url]
        self.assertEqual(len(asked), 1)
        self.assertIn("1TUP", asked[0])
        cmd.delete("*")
        self.assertEqual(
            cmd.fetch_uniprot_ligands("P04637", 0, 1, "lig"), ["P04637_AZM"]
        )
        self.assertIn("P04637_ZN", cmd.fetch_uniprot_ligands("P04637", 0, 0, "any"))
        self.assertRaises(CmdException, cmd.fetch_uniprot_ligands, "P04637", 1000)

    def test_fold_alphafold(self):
        entries = [
            {"entryId": "AF-P04637-F1", "pdbUrl": "http://x/F1.pdb"},
            {"entryId": "AF-P04637-F2", "pdbUrl": "http://x/F2.pdb"},
        ]
        self.serve(
            {
                "prediction/P04637": json.dumps(entries),
                "F1.pdb": model([30, 95]),
                "F2.pdb": model([60]),
            }
        )
        self.assertEqual(cmd.fetch_fold("p04637"), "AF_P04637_F1")
        self.assertEqual(cmd.count_atoms("AF_P04637_F1 and name CA"), 2)
        self.assertEqual(len(colors("AF_P04637_F1")), 2)
        self.assertIn("plddt_legend", cmd.get_names())
        self.assertEqual(cmd.fetch_fold("AF-P04637-F2", "second", 0), "second")
        self.assertEqual(cmd.count_atoms("second and name CA"), 1)
        self.assertRaises(CmdException, cmd.fetch_fold, "not a protein!")

    def test_fold_esmfold(self):
        # an accession AlphaFold DB does not have: sequence from UniProt
        self.serve(
            {
                "prediction/Q9XXX1": "",
                "Q9XXX1.fasta": ">sp|Q9XXX1|X\nACD\nEF\n",
                "esmatlas": model([0.3, 0.6, 0.8, 0.95, 0.5]),
            }
        )
        self.assertEqual(cmd.fetch_fold("Q9XXX1"), "ESM_Q9XXX1")
        self.assertEqual(self.asked[-1][1], b"ACDEF")
        # the server's 0 to 1 becomes 0 to 100, one color per bin
        bfactors = []
        cmd.iterate("ESM_Q9XXX1 and name CA", "bfactors.append(b)", space=locals())
        self.assertEqual([round(b) for b in bfactors], [30, 60, 80, 95, 50])
        self.assertEqual(len(colors("ESM_Q9XXX1")), 4)
        name = cmd.fetch_fold("acd ef", color=0)  # a sequence, spaces allowed
        self.assertTrue(name.startswith("ESM_"))
        self.assertEqual(self.asked[-1][1], b"ACDEF")
        self.serve({"esmatlas": "Internal Server Error"})
        self.assertRaises(CmdException, cmd.fetch_fold, "ACDEF")

    def test_color_plddt(self):
        cmd.fab("AAAA", "m")
        for resi, b in enumerate((49.9, 50, 70, 90), 1):
            cmd.alter("resi %d" % resi, "b=%r" % b)
        cmd.color_plddt("m", legend=0)
        self.assertNotIn("plddt_legend", cmd.get_names())
        for resi, name in enumerate(row[0] for row in fetching.PLDDT):
            self.assertEqual(
                colors("m and resi %d" % (resi + 1)),
                {cmd.get_color_tuple(name)},
                name,
            )
        cmd.color_plddt()
        cmd.color_plddt()  # the legend is replaced, not doubled
        self.assertEqual(cmd.count_atoms("plddt_legend"), 4)
        self.assertRaises(CmdException, cmd.color_plddt, "none")

    def test_registered(self):
        for name in (
            "fetch_by_uniprot",
            "fetch_fold",
            "fetch_uniprot_ligands",
            "list_uniprot_structures",
            "color_plddt",
        ):
            self.assertIn(name, cmd.keyword)
            self.assertIn("DESCRIPTION", getattr(cmd, name).__doc__)

    @testing.requires("network", "no_run_all")
    def test_live(self):
        self.assertIn("1TUP", cmd.list_uniprot_structures("P04637"))
        self.assertEqual(len(cmd.fetch_uniprot_ligands("P00918", limit=3)), 3)
        self.assertEqual(cmd.fetch_fold("AF-P04637-F1"), "AF_P04637_F1")
        self.assertGreater(cmd.count_atoms("AF_P04637_F1"), 1000)
        name = cmd.fetch_fold("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
        bfactors = []
        cmd.iterate(name, "bfactors.append(b)", space=locals())
        self.assertGreater(max(bfactors), 1.5)  # 0 to 100, not 0 to 1
