# -c

/print "BEGIN-LOG"

load dat/odd01.pdb
iterate all, print name

save cmp/C0700odd.1.pdb

set retain_order
save cmp/C0700odd.2.pdb

set pdb_no_end_record
save cmp/C0700odd.3.pdb

set pdb_retain_ids
save cmp/C0700odd.4.pdb

iterate all,print name

reinit

set pdb_literal_names
set retain_order
load dat/odd01.pdb
iterate all,print name
save cmp/C0700odd.5.pdb

dele all

load cmp/C0700odd.5.pdb
save cmp/C0700odd.6.pdb

reinit

load cmp/C0700odd.5.pdb
iterate all,print name
save cmp/C0700odd.7.pdb

/print "END-LOG"