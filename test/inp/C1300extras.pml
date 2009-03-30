# -c

/print "BEGIN-LOG"
delete all

fragment methane
remove hydro and not first hydro
alter all,segi='METH'
alter hydro, alt='A'
alter hydro, text_type='Htype'
alter all, ss='L'
alter all, b=rank
alter all, q=index
alter all, elec_radius = vdw+1.0

# alternate label2 routine (ugly and quite limited, but bypasses python eval)

cmd.label2("all","'some text'")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","'some text '+'more text'")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all"," 'some text ' + 'more text' ")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","model")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","resv")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","model+' index:'+index+ ' type:' +type")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","segi+' '+chain+' '+resi+'`'+resn+' '+str(resv)+name+'`'+alt")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","str(index)+' '+str(rank)+' '+str(ID)")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","flags")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","str(numeric_type)+text_type")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","str(formal_charge) + ' ' + str(partial_charge)")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all","ss + str(cartoon) + str(color)")
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all",'str(valence) + " " + str(geom)')
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

cmd.label2("all",'str(elec_radius)+" "+str(vdw)+" "+b+" "+q+" "+elem')
cmd.iterate("all",r"print name +' -> ' + '\'' + label + '\''")

print "END-LOG"