"""
Find a sequence motif in proteins. "pepseq" selects an exact sequence;
this takes a regular expression. New code, MIT.
"""

import re

import pymol

cmd = __import__("sys").modules["pymol.cmd"]


def search_sequence(pattern, name="", selection="all", quiet=1, *, _self=cmd):
    """
    DESCRIPTION

        "search_sequence" selects the residues of every protein chain whose
        sequence matches a pattern. The pattern is a sequence in one-letter
        code or a regular expression.

    USAGE

        search_sequence pattern [, name [, selection ]]

    ARGUMENTS

        pattern = string: for example GDSL, or N[^P][ST] for a
        glycosylation site ("." is any residue)

        name = string: name of the selection {default: seq_ and the letters
        of the pattern}

        selection = string: where to look {default: all}

    EXAMPLE

        search_sequence HPGG
        search_sequence N[^P][ST], sequon, chain A

    NOTES

        Returns the number of matches. A match does not continue from one
        chain into the next. Upper and lower case are the same. For an
        exact sequence, "select name, pepseq GDSL" does the same.

    SEE ALSO

        select, wizard find
    """
    try:
        regex = re.compile(pattern, re.IGNORECASE)
    except re.error as error:
        raise pymol.CmdException("pattern '%s': %s" % (pattern, error)) from None
    chains = {}  # (object, segi, chain) -> [(atom index, one-letter code)]
    _self.iterate(
        "(%s) and polymer.protein and name CA and alt +A" % selection,
        "chains.setdefault((model, segi, chain), []).append((index, oneletter))",
        space={"chains": chains},
    )
    name = name or "seq_" + re.sub(r"\W", "", pattern)
    _self.select(name, "none")
    count = 0
    for (model, _, _), residues in chains.items():
        sequence = "".join(letter or "X" for _, letter in residues)
        found = []
        for match in regex.finditer(sequence):
            found += [index for index, _ in residues[match.start() : match.end()]]
            count += match.end() > match.start()
        if found:
            _self.select_list("_seq_hit", model, found, mode="index")
            _self.select(name, "?%s or byres _seq_hit" % name)
            _self.delete("_seq_hit")
    if not int(quiet):
        print(
            " search_sequence: %d matches, %d residues in '%s'"
            % (count, _self.count_atoms("?%s and name CA" % name), name)
        )
    return count
