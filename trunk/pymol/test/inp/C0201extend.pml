# -c

/print "BEGIN-LOG"

def foo(moo=2): print moo
cmd.extend('foo',foo)

foo
foo 3
foo moo=5
foo ?

def fie(moo=3):\
   "help for fie" \
   print moo

print fie.__doc__

help fie
cmd.extend('fie',fie)
help fie

/print "END-LOG"
