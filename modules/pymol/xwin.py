
def machine_get_clipboard():
   result = []
   from popen2 import popen2
   pipe = popen2("/usr/bin/wish")
   pipe[1].write("puts [ selection get ]\n")
   pipe[1].write("exit\n")
   pipe[1].close()
   while 1:
      l=pipe[0].readline()
      if not l: break
      result.append(l)
   return result

                 
