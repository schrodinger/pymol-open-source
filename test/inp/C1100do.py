# -c

from pymol import cmd

print "BEGIN-LOG"

cmd.do("""
load dat/pept.pdb
show sticks,pept
set ray_default_renderer=2
ray
dele all
""")

cmd.do(
   ['load dat/il2.pdb',
    'show dots,il2',
    'ray',
    'dele all'
    ])


cmd.do('load dat/3al1.pdb\nshow sph,3al1\nray\ndele all')

cmd.do('fragment arg')
cmd.do('ray')

# see if when can embed Python blocks inside of PyMOL inside of Python

cmd.do('for a in range(1,10):\\\n   print a\\\n   print a+10')

cmd.do(r"""
for a in range(11,21):\
  print a\
  print a+10
  """)


print "END-LOG"





