from chempy.champ import Champ

ch = Champ()
p1 = ch.insert_pattern_string("C(N)(F)O")
ch.pattern_dump(p1)


ch = Champ()
p1 = ch.insert_pattern_string("C1(N)(F).O1")
ch.pattern_dump(p1)

