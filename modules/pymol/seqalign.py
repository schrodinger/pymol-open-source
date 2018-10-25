'''
Sequence alignment stuff.

Based on psico.seqalign

(c) 2018 Schrodinger, Inc.
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def needle_alignment(s1, s2):
    '''
DESCRIPTION

    Does a Needleman-Wunsch Alignment of sequence s1 and s2 and
    returns a Bio.Align.MultipleSeqAlignment object.
    '''
    from Bio import pairwise2
    from Bio.Align import MultipleSeqAlignment
    from Bio.SubsMat.MatrixInfo import blosum62

    def match_callback(c1, c2):
        return blosum62.get((c1, c2), 1 if c1 == c2 else -4)

    alns = pairwise2.align.globalcs(s1, s2,
            match_callback, -10., -.5,
            one_alignment_only=True)

    a = MultipleSeqAlignment([])
    a.add_sequence("s1", alns[0][0])
    a.add_sequence("s2", alns[0][1])
    return a


def alignment_mapping(seq1, seq2):
    '''
DESCRIPTION

    Returns an iterator with seq1 indices mapped to seq2 indices

    >>> mapping = dict(alignment_mapping(s1, s2))
    '''
    i, j = -1, -1
    for a, b in zip(seq1, seq2):
        if a != '-': i += 1
        if b != '-': j += 1
        if a != '-' and b != '-': yield i, j


def aln_magic_format(infile):
    '''
DESCRIPTION

    Guess alignment file format.
    '''
    with open(infile) as handle:
        for line in handle:
            if len(line.rstrip()) > 0:
                break
    if line.startswith('CLUSTAL') or line.startswith('MUSCLE'):
        informat = 'clustal'
    elif line.startswith('>P1;'):
        informat = 'pir'
    elif line.startswith('>'):
        informat = 'fasta'
    elif line.startswith('# STOCKHOLM'):
        informat = 'stockholm'
    elif line.startswith('Align '):
        informat = 'fatcat'
    elif line.startswith('# ProSMART Alignment File'):
        informat = 'prosmart'
    else:
        informat = 'emboss'
    return informat


def aln_magic_read(infile, format=None):
    '''
DESCRIPTION

    Wrapper for Bio.AlignIO.read that guesses alignment file format.
    '''
    from Bio import AlignIO

    if not format:
        format = aln_magic_format(infile)

    with open(infile) as handle:
        return AlignIO.read(handle, format)


def load_aln_multi(filename, object=None, mapping='', alignioformat='', guide='',
        quiet=1, _self=cmd):
    '''
DESCRIPTION

    Load a multiple sequence alignment from file and apply it to already
    loaded structures.

    Requires Biopython (conda install biopython)

USAGE

    load_aln_multi filename [, object [, mapping [, alignioformat [, guide ]]]]

ARGUMENTS

    filename = str: alignment file

    object = str: name of the new alignment object {default: filename prefix}

    mapping = str: space separated list of sequence-id to object-name
    mappings {default: assume sequence-id = object-name}

    alignioformat = str: file format, see http://biopython.org/wiki/AlignIO
    {default: guess from first line in file}

    guide = str: object name of "guide" object (PyMOL only knows many-to-one
    alignment, not true multiple alignments) {default: random object}

SEE ALSO

    psico.importing.load_aln
    '''
    import os

    quiet = int(quiet)
    if object is None:
        object = os.path.basename(filename).rsplit('.', 1)[0]

    # load alignment file
    alignment = aln_magic_read(cmd.exp_path(filename), alignioformat)

    # object name mapping
    if isinstance(mapping, str):
        mapping = mapping.split()
        mapping = dict(zip(mapping[0::2], mapping[1::2]))

    i2index = [None] * len(alignment)
    onames = [None] * len(alignment)
    indices = [-1] * len(alignment)

    # get structure models and sequences
    for r, record in enumerate(alignment):
        oname = mapping.get(record.id, record.id)

        if not oname:
            continue

        atoms = []
        n = _self.iterate('?{} & guide'.format(oname), 'atoms.append((oneletter or "?", index))',
                space=locals())

        if n == 0:
            print(" Warning: no atoms for object '{}'".format(oname))
            continue

        # align sequences from file to structures
        aln = needle_alignment(str(record.seq), ''.join(a[0] for a in atoms))

        # get index mappings
        i2index[r] = dict((i, atoms[j][1])
                for (i, j) in alignment_mapping(*aln))
        onames[r] = oname

    if len(onames) < 2:
        raise CmdException('Failed to map alignment to objects')

    # build alignment list
    raw = []
    for c in range(len(alignment[0])):
        rawcol = []
        for r, aa in enumerate(alignment[:, c]):
            if onames[r] is None:
                continue
            if aa != '-':
                indices[r] += 1
                index = i2index[r].get(indices[r], -1)
                if index >= 0:
                    rawcol.append((onames[r], index))
        if len(rawcol) > 1:
            raw.append(rawcol)

    _self.set_raw_alignment(object, raw, guide=guide)
    return raw


# vi:expandtab:smarttab
