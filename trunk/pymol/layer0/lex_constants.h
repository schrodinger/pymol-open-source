/*
 * This file contains source code for the PyMOL computer program
 * Copyright (c) Schrodinger, LLC.
 *
 * Registration of all string "globals" in PyMOLGlobals::Lexicon
 *
 * This file must only be included at its two designated places.
 */

#ifdef LEX_CONSTANTS_IMPL
#undef LEX_CONSTANTS_IMPL
#define LEX_CONST(x) G->lex_const.x = LexIdx(G, #x);
#else
#define LEX_CONST(x) x,
#endif

LEX_CONST(ARG)
LEX_CONST(ASN)
LEX_CONST(ASP)
LEX_CONST(C)
LEX_CONST(CA)
LEX_CONST(CB)
LEX_CONST(CD)
LEX_CONST(GLN)
LEX_CONST(GLU)
LEX_CONST(HOH)
LEX_CONST(LEU)
LEX_CONST(LYS)
LEX_CONST(MET)
LEX_CONST(MSE)
LEX_CONST(N)
LEX_CONST(O)
LEX_CONST(OXT)
LEX_CONST(P)
LEX_CONST(PHE)
LEX_CONST(PRO)
LEX_CONST(TYR)
LEX_CONST(H1)
LEX_CONST(H3)

#undef LEX_CONST
