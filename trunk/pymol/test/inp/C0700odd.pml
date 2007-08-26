# -c

/print "BEGIN-LOG"

# This file tests PyMOL's ability to handle odd PDB files, 
# both in terms of reading and writing intact, but also converting to PDB compliant format
# as well as converting to AMBER atom naming format

# it also test PyMOL's ability to retain identifiers and ordering

unset pdb_hetatm_guess_valences
reinitialize store_defaults

# test default preservation 

load dat/odd01.pdb
iterate all, print name
save cmp/C0700odd.1.pdb

# test order retention after read

set retain_order
iterate all, print name
save cmp/C0700odd.2.pdb

# test END record control

set pdb_no_end_record
save cmp/C0700odd.3.pdb

# test retention of IDs

set pdb_retain_ids
save cmp/C0700odd.4.pdb

reinit

# test order retention before read

set retain_order
load dat/odd01.pdb
iterate all,print name

# test convert AMBER->PDB on write

set pdb_reformat_names_mode, 1
save cmp/C0700odd.5.pdb

dele all

# test convert AMBER->PDB on read

load cmp/C0700odd.3.pdb

set pdb_reformat_names_mode, 0
save cmp/C0700odd.6.pdb
iterate all,print name

reinit

# test AMBER->PDB on both

set pdb_reformat_names_mode, 1

load cmp/C0700odd.3.pdb
iterate all,print name
save cmp/C0700odd.7.pdb

# test PDB->AMBER conversion on write

set pdb_reformat_names_mode,2
save cmp/C0700odd.8.pdb

dele all

# test PDB->AMBER convert on read
load cmp/C0700odd.5.pdb
iterate all,print name

set pdb_reformat_names_mode,0
save cmp/C0700odd.9.pdb

# test AMBER->PDB back conversion on write

set pdb_reformat_names_mode,1
save cmp/C0700odd.A.pdb

# test literal names override on write 

set pdb_literal_names, 1
save cmp/C0700odd.B.pdb
dele all

# test literal names override on read

load cmp/C0700odd.3.pdb
iterate all,print name
dele all
load cmp/C0700odd.5.pdb
iterate all,print name

# check halogen handling

reinit

load dat/small02.pdb
save cmp/C0700odd.C.pdb

reinit

set pdb_reformat_names_mode, 1
load dat/small02.pdb
iterate all, print name
save cmp/C0700odd.D.pdb

set pdb_reformat_names_mode, 2
load dat/small02.pdb
iterate all, print name
save cmp/C0700odd.E.pdb

reinit

# test literal names preservation capability

set pdb_literal_names
set retain_order
set pdb_retain_ids
set pdb_no_end_record
unset pdb_use_ter_records

load dat/odd01.pdb
save cmp/C0700odd.F.pdb

reinit

load dat/helix_amber.pdb
iterate resi 2, print name, elem

reinit
load dat/odd02.pdb
iterate all, print name, elem
save cmp/C0700.odd.G.pdb
set pdb_reformat_names_mode,1
save cmp/C0700.odd.H.pdb

/print "END-LOG"