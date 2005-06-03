

import os

def machine_get_clipboard():
    result = []
    from popen2 import popen2
    wish_path = os.environ['PYMOL_PATH']+"/ext/bin/wish8.0" # dubious...
    if not os.path.exists(wish_path):
        wish_path = "/usr/bin/wish"
    if os.path.exists(wish_path):
        pipe = popen2(wish_path)
        pipe[1].write("puts [ selection get ]\n")
        pipe[1].write("exit\n")
        pipe[1].close()
        while 1:
            l=pipe[0].readline()
            if not l: break
            result.append(l)
    return result

                      
    
