# -c

set auto_number_selections
/print "BEGIN-LOG"

load dat/pept.pdb

select (resi trp)
indicate (name ca)
deselect
select 1/o
select 1/o+ca+n
select /pept

select bb = (name ca or name c or name n)
select bb = (name ca,c,n)
select bb = */ca+c+n
select bb = (*/ca,c,n)
select bb = (not (not (name ca or name c or name n)))
select bb = (not (not (n. c|n. ca|n. n)))
select bb = (not (not (not (not (not (not (not (not (not (not (not (not (not (not (*/ca,c,n)))))))))))))))
select bb = (all and all and all and (not (not all)) and (not none) and all and not (not (n;ca,c,n)))
select bb,(name ca or (n;c | (n. n) | (*/n)))
select bb, \
name ca \
or name c \
or name n

select bb = bb
select bb = not not bb
select cc = not bb
select bb = not cc

enable
enable bb
enable cc
disable cc
hide lines,cc
hide lines,bb
show (bb)
set ray_default_renderer=2
ray

dele all

load dat/1tii.pdb
count_atoms 1+2+3+4+5+6+7+8+9+10+11+12+13+14+15+16+17+18+19+20+21+22+23+24+25+26+27+28+29+30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57+58+59+60+61+62+63+64+65+66+67+68+69+70+71+72+73+74+75+76+77+78+79+80+81+82+83+84+85+86+87+88+89+90+91+92+93+94+95+96+97+98+99+100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118+119+120+121+122+123+124+125+126+127+128+129+130+131+132+133+134+135+136+137+138+139+140+141+142+143+144+145+146+147+148+149+150+151+152+153+154+155+156+157+158+159+160+161+162+163+164+165+166+167+168+169+170+171+172+173+174+175+176+177+178+179+180+181+182+183+184+185+186+187+195+196+197+198+199+200+201+202+203+204+205+206+207+208+209+210+211+212+213+214+215+216+217+218+219+220+221+222+223+224+225+226+227+228+229+230/

print cmd.count_atoms("resi 1+2+3+4+5+6+7+8+9+10+11+12+13+14+15+16+17+18+19+20+21+22+23+24+25+26+27+28+29+30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57+58+59+60+61+62+63+64+65+66+67+68+69+70+71+72+73+74+75+76+77+78+79+80+81+82+83+84+85+86+87+88+89+90+91+92+93+94+95+96+97+98+99+100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118+119+120+121+122+123+124+125+126+127+128+129+130+131+132+133+134+135+136+137+138+139+140+141+142+143+144+145+146+147+148+149+150+151+152+153+154+155+156+157+158+159+160+161+162+163+164+165+166+167+168+169+170+171+172+173+174+175+176+177+178+179+180+181+182+183+184+185+186+187+195+196+197+198+199+200+201+202+203+204+205+206+207+208+209+210+211+212+213+214+215+216+217+218+219+220+221+222+223+224+225+226+227+228+229+230")

# confirm ability to support 'model' as an object name sequence

dele all
load dat/pept.pdb
load dat/tiny.pdb, model
load dat/small01.mol, test+me
load dat/small02.pdb, mymodel
load dat/small03.mol2, model2

# should error
count_atoms model

# shouldn't error
count_atoms model model
count_atoms /model
count_atoms pept
count_atoms model model pept
count_atoms /model+pept
count_atoms /pept+model
count_atoms /test+\me
count_atoms /pept+model+test\+me
count_atoms test+me
count_atoms test+me model model
count_atoms mymodel
count_atoms model2

/print "END-LOG"









