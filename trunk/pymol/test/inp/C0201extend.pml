# -c

/print "BEGIN-LOG"

def foo(moo=2): print moo
cmd.extend('foo',foo)

foo
foo 3
foo moo=5
foo ?

/print "END-LOG"
