
aa_dict = {
   'A' : 'ala',
   'C' : 'cys',
   'D' : 'asp',
   'E' : 'glu',
   'F' : 'phe',

   'G' : 'gly',
   'H' : 'his',
   'I' : 'ile',
   'K' : 'lys',
   'L' : 'leu',

   'M' : 'met',
   'N' : 'asn',
   'P' : 'pro',
   'Q' : 'gln',
   'R' : 'arg',

   'S' : 'ser',
   'T' : 'thr',
   'V' : 'val',
   'W' : 'trp',
   'Y' : 'tyr',
}

from pymol import editor
from pymol import cmd

def build(object_name, sequence, first_residue = "1"):
   if len(sequence):
      code = sequence[0]
      cmd.fragment(aa_dict[code],object_name)
      cmd.alter(object_name,'resi="%s"'%first_residue)
      cmd.edit(object_name+" and name C")
      for code in sequence[1:]:
         editor.attach_amino_acid("pk1",aa_dict[code])
      cmd.edit()
      
build("poly_ala","ACDEFGHIKLMNPQRSTVWY")

cmd.zoom()

   
   

