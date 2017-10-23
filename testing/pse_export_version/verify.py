from pymol import cmd

def myasserttrue(b):
    if not b:
        raise UserWarning('myasserttrue')

def myassertequal(a, b):
    if a != b:
        raise UserWarning('myassertequal %s != %s' % (a, b))

myassertequal(cmd.get_names(), ['1rx1', '1rx1_2fofc', 'mesh'])
myassertequal(cmd.count_atoms('color forest'), 1268)
myassertequal(cmd.count_atoms('color red'), 48)
myassertequal(cmd.count_atoms('color orange'), 1)
myasserttrue(cmd.get('bg_rgb') in ('yellow', '[ 1.00000, 1.00000, 0.00000 ]', '0xffff00'))
myassertequal(cmd.get_extent('mesh'), [[1.3251924514770508, 15.123332977294922, -12.337624549865723], [51.682502746582031, 73.096115112304688, 37.698299407958984]])
myassertequal(cmd.get_wizard_stack()[0].__class__.__name__, 'Mutagenesis')

cmd.scene('F1', 'recall')
myassertequal(cmd.get_wizard_stack()[-1].message[0], "First Scene")
