

import os

def machine_get_clipboard():
    result = []
    wish_path = os.environ['PYMOL_PATH']+"/ext/bin/wish8.0" # dubious...
    if not os.path.exists(wish_path):
        wish_path = "/usr/bin/wish"
    if os.path.exists(wish_path):
        import subprocess
        p = subprocess.Popen([wish_path],
                universal_newlines=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
        p.stdin.write("puts [ selection get ]\n")
        p.stdin.write("exit\n")
        p.stdin.close()
        result = p.stdout.readlines()
    return result
