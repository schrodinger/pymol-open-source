"""
Structures by UniProt accession: the entries of the PDB, predicted models
from AlphaFold DB or ESMFold, and pLDDT colors. Ported from the pymolish
plugin (MIT, same author); new code, MIT. Standard library HTTP only.
"""

import hashlib
import json
import re
import urllib.error
import urllib.request

import pymol

cmd = __import__("sys").modules["pymol.cmd"]

UNIPROT = "https://rest.uniprot.org/uniprotkb/"
ALPHAFOLD = "https://alphafold.ebi.ac.uk/api/prediction/"
ESMFOLD = "https://api.esmatlas.com/foldSequence/v1/pdb/"

_ACCESSION = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)
_ALPHAFOLD_ID = re.compile(r"AF-([A-Z0-9]+)-F(\d+)")

# AlphaFold DB confidence colors: name, selection, legend label, RGB
PLDDT = (
    ("plddt_very_low", "b < 50", "<50", (1.00, 0.49, 0.27)),
    ("plddt_low", "b >= 50 and b < 70", "50-70", (1.00, 0.86, 0.07)),
    ("plddt_confident", "b >= 70 and b < 90", "70-90", (0.40, 0.80, 0.95)),
    ("plddt_very_high", "b >= 90", "90-100", (0.00, 0.33, 0.84)),
)
LEGEND = "plddt_legend"


def _get(url, data=None):
    """Text of the answer; "" when the server does not know the entry."""
    request = urllib.request.Request(
        url, data=data, headers={"User-Agent": "pymolish-studio"}
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as answer:
            return answer.read().decode()
    except urllib.error.HTTPError as error:
        if error.code in (400, 404):
            return ""
        raise pymol.CmdException("%s: %s" % (url, error)) from None
    except OSError as error:
        raise pymol.CmdException("%s: %s" % (url, error)) from None


def pdb_entries(accession):
    """(PDB id, method, resolution, chains) of every PDB entry that UniProt
    lists for an accession."""
    accession = accession.strip().upper()
    text = _get(UNIPROT + accession + ".json?fields=xref_pdb")
    if not text:
        raise pymol.CmdException("unknown UniProt accession '%s'" % accession)
    rows = []
    for ref in json.loads(text).get("uniProtKBCrossReferences", []):
        if ref.get("database") == "PDB":
            about = {p["key"]: p["value"] for p in ref.get("properties", [])}
            rows.append(
                (
                    ref["id"],
                    about.get("Method", ""),
                    about.get("Resolution", ""),
                    about.get("Chains", ""),
                )
            )
    return rows


def list_uniprot_structures(accession, quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "list_uniprot_structures" prints the PDB entries of a protein, with
        method, resolution and the chains that cover it. Nothing is loaded.

    USAGE

        list_uniprot_structures accession

    EXAMPLE

        list_uniprot_structures P04637

    SEE ALSO

        fetch_by_uniprot, fetch_fold
    """
    rows = pdb_entries(accession)
    if not int(quiet):
        for row in rows:
            print(" %-5s %-8s %-8s %s" % row)
        print(" list_uniprot_structures: %d entries" % len(rows))
    return [row[0] for row in rows]


def fetch_by_uniprot(accession, limit=0, group="", async_=0, quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "fetch_by_uniprot" downloads the PDB entries of a protein, given
        its UniProt accession, into one group.

    USAGE

        fetch_by_uniprot accession [, limit [, group [, async_ ]]]

    ARGUMENTS

        accession = string: UniProt accession, for example P04637

        limit = integer: load only the first entries {default: 0, all}

        group = string: name of the group {default: the accession}

        async_ = 0/1: download in the background {default: 0}

    EXAMPLE

        fetch_by_uniprot P00918, 5

    NOTES

        The objects are named <accession>_<PDB id>. The files go where
        "fetch" puts them, see the "fetch_path" setting. A well studied
        protein has hundreds of entries: look at "list_uniprot_structures"
        first, or give a limit.

    SEE ALSO

        list_uniprot_structures, fetch_fold, fetch
    """
    if int(async_):
        return _self.async_(
            fetch_by_uniprot, accession, limit, group, 0, quiet, _self=_self
        )
    accession = accession.strip().upper()
    codes = [row[0] for row in pdb_entries(accession)]
    if not codes:
        raise pymol.CmdException("UniProt lists no PDB entry for " + accession)
    codes = codes[: int(limit) or None]
    names = []
    for code in codes:
        name = "%s_%s" % (accession, code)
        try:
            _self.fetch(code, name, zoom=0, quiet=1)
        except pymol.CmdException as error:
            print(" fetch_by_uniprot: skipped %s: %s" % (code, error))
            continue
        names.append(name)
    if not names:
        raise pymol.CmdException("no entry could be downloaded")
    _self.group(group or accession, " ".join(names))
    _self.zoom(group or accession)
    if not int(quiet):
        print(" fetch_by_uniprot: loaded %d of %d" % (len(names), len(codes)))
    return names


def fetch_fold(query, name="", color=1, async_=0, quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "fetch_fold" loads a predicted structure. A UniProt accession or an
        AlphaFold id is looked up in AlphaFold DB; a protein sequence is
        folded by the public ESMFold server. An accession that AlphaFold DB
        does not have is folded by ESMFold from its UniProt sequence.

    USAGE

        fetch_fold query [, name [, color [, async_ ]]]

    ARGUMENTS

        query = string: UniProt accession (P04637), AlphaFold id
        (AF-P04637-F1) or sequence in one-letter code

        name = string: object name {default: AF_<accession>_F<n> or
        ESM_<accession or hash>}

        color = 0/1: color by confidence, see "color_plddt" {default: 1}

        async_ = 0/1: download in the background {default: 0}

    EXAMPLE

        fetch_fold P04637
        fetch_fold MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ, peptide

    NOTES

        The sequence is sent to api.esmatlas.com. The server takes about
        400 residues at most. The confidence (pLDDT, 0 to 100) is in the
        B-factor column.

    SEE ALSO

        color_plddt, fetch_by_uniprot
    """
    if int(async_):
        return _self.async_(fetch_fold, query, name, color, 0, quiet, _self=_self)
    query = "".join(query.split()).upper()
    known = _ALPHAFOLD_ID.match(query)
    accession, fragment = (known.group(1), known.group(2)) if known else (query, "1")
    pdb = sequence = ""
    if _ACCESSION.fullmatch(accession):
        models = json.loads(_get(ALPHAFOLD + accession) or "[]")
        wanted = "AF-%s-F%s" % (accession, fragment)
        if models:
            model = next((m for m in models if m.get("entryId") == wanted), models[0])
            pdb = _get(model["pdbUrl"])
            name = name or model["entryId"].replace("-", "_")
        else:
            fasta = _get(UNIPROT + accession + ".fasta")
            sequence = "".join(fasta.splitlines()[1:])
            name = name or "ESM_" + accession
    elif query and set(query) <= set("ACDEFGHIKLMNPQRSTVWY"):
        sequence = query
        name = name or "ESM_" + hashlib.sha1(query.encode()).hexdigest()[:6]
    else:
        raise pymol.CmdException(
            "give a UniProt accession, an AlphaFold id or a protein sequence"
        )
    if sequence:
        pdb = _get(ESMFOLD, sequence.encode())
    if "ATOM" not in pdb:
        raise pymol.CmdException("no model for '%s'" % query[:40])
    _self.read_pdbstr(pdb, name)
    if sequence:  # the ESMFold server writes pLDDT as 0 to 1
        _self.alter("%" + name, "b = b * 100")
    if int(color):
        color_plddt("%" + name, _self=_self)
    if not int(quiet):
        print(" fetch_fold: loaded '%s'" % name)
    return name


def color_plddt(selection="all", legend=1, quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "color_plddt" colors a predicted model by confidence with the four
        AlphaFold DB colors: orange below 50, yellow to 70, light blue to
        90, blue above. The confidence (pLDDT) is read from the B-factor
        column, where AlphaFold DB and "fetch_fold" put it.

    USAGE

        color_plddt [ selection [, legend ]]

    ARGUMENTS

        selection = string: atoms to color {default: all}

        legend = 0/1: put four labeled spheres next to the model, object
        "plddt_legend" {default: 1}

    NOTES

        For a continuous scale use
        "spectrum b, red_yellow_blue, minimum=0, maximum=100".

    SEE ALSO

        fetch_fold, spectrum
    """
    if not _self.count_atoms(selection):
        raise pymol.CmdException("no atoms in '%s'" % selection)
    for name, test, _, rgb in PLDDT:
        _self.set_color(name, rgb)
        _self.color(name, "(%s) and %s" % (selection, test))
    _self.delete(LEGEND)
    if int(legend):
        (_, _, zmin), (xmax, ymax, _) = _self.get_extent(selection)
        for row, (name, _, label, _) in enumerate(PLDDT):
            _self.pseudoatom(
                LEGEND,
                pos=[xmax + 10, ymax - 5 * row, zmin],
                label=label,
                color=name,
            )
        _self.show_as("spheres", LEGEND)
        _self.show("labels", LEGEND)
        _self.set("sphere_scale", 1.5, LEGEND)
        _self.set("label_position", [3.5, 0, 0], LEGEND)
    if not int(quiet):
        print(" color_plddt: colored %d atoms" % _self.count_atoms(selection))
