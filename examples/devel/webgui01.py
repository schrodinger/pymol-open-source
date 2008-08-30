import SocketServer
import BaseHTTPServer

import time
import cgi
import threading

from pymol import cmd

_server = None
def _shutdown(self_cmd=cmd):
    if _server != None:
        _server.socket.close()
    self_cmd.quit()

# Note, this handler assumes that PyMOL is running as a global singleton

def get_status(out, self_cmd=cmd):
    out.write('<html><body>\n')
    out.write('<h3>PyMOL WebGUI Proof of Concept</h3>\n')
    out.write('<table><tr>\n')
    out.write('<td><form action="./status.pymol"><button type="submit">Refresh</button></form></td>\n')
    out.write('<td><form target="_new" action="./ray.pymol"><button type="submit">Ray</button></form></td>\n')    
    out.write('<td><form action="./quit.pymol"><button type="submit">Quit</button></form></td>\n')
    out.write('</tr></table>')
    out.write('<a href="./status.pymol?load">load $TUT/1hpv.pdb</a>\n')
    names = self_cmd.get_names('objects')
    if not len(names):
        out.write('<p>No objects loaded.</p>\n')
    else:
        out.write('<p>Loaded Objects:</p><ul>\n')
        for name in names:
            out.write('<li>%s</li>\n'%name)
        out.write('</ul>\n')
    out.write('</body></html>\n')

def get_start(out, self_cmd=cmd):
    window_open="javascript: window.open('./status.pymol','hello', 'location=no,status=no,toolbar=no,width=400,height=600,top=0,left=880');"
    out.write('<html><body onload="'+window_open+'">\n')
    out.write('<a href="./start.pymol" onclick="'+window_open+'">Launch WebGUI\n')
    out.write('</body></html>')

def get_ray(out, self_cmd=cmd):
    
    self_cmd.ray()
    self_cmd.png('tmp.png')
    self_cmd.sync()
    # need to work on synchronization...
    out.write(open('tmp.png').read())
    
class PymolHandler(BaseHTTPServer.BaseHTTPRequestHandler):

    def do_pymol(self):
        if "ray.pymol" in self.path: # send image
            self.send_response(200)
            self.send_header('Content-type',	'image/x-png')
            self.end_headers()
            get_ray(self.wfile)
        else:
            if "load" in self.path: # load a structure
                cmd.load("$TUT/1hpv.pdb")
            self.send_response(200)
            self.send_header('Content-type',	'text/html')
            self.end_headers()
            if "status.pymol" in self.path:
                get_status(self.wfile)
            elif "quit.pymol" in self.path:
                self.wfile.write('<html><body><p>Quitting...</p></body></html>')
                self.wfile.flush()
                _shutdown()
            else: # start page
                get_start(self.wfile)
        self.wfile.flush()
        
    def do_GET(self):
        try:
            doc = self.path.split('?')[0]
            if doc.endswith('.pymol'): # PyMOL
                self.do_pymol()
            elif doc.endswith('.html'):
                f = open('.'+self.path) # UNSAFE!!!
                self.send_response(200)
                self.send_header('Content-type',	'text/html')
                self.end_headers()
                self.wfile.write(f.read())
                f.close()
        except IOError:
            self.send_error(404,'File Not Found: %s' % self.path)

    def do_POST(self):
        global rootnode
        try:
            ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
            if ctype == 'multipart/form-data':
                query=cgi.parse_multipart(self.rfile, pdict)
            self.send_response(301)
            
            self.end_headers()
            upfilecontent = query.get('upfile')
            print "filecontent", upfilecontent[0]
            self.wfile.write('<HTML>POST OK.<BR><BR>');
            self.wfile.write(upfilecontent[0]);
            
        except :
            pass

class ThreadingHTTPServer(SocketServer.ThreadingMixIn,
                          BaseHTTPServer.HTTPServer):
    pass

def main():
    try:
        global _server
        _server = ThreadingHTTPServer(('', 8080), PymolHandler)
        print 'started httpserver...'
        _server.serve_forever()
    except KeyboardInterrupt:
        print '^C received, shutting down server'
        _server.socket.close()

def open_browser():
    time.sleep(1)
    import os
    os.system('open http://localhost:8080/start.pymol')

if __name__ == '__main__':
    main()

if __name__ == 'pymol':
    t = threading.Thread(target=main)
    t.setDaemon(1)
    t.start()

    t = threading.Thread(target=open_browser)
    t.setDaemon(1)
    t.start()

    
"""

okay, what we need now is a simple safe way to pass comands and
arguments through the client brower, with escape characters, etc.

ideally, we'd like a javascript object which responds to the same messages as PyMOL

for our initial test, let's just put up some button links to control representations.

also outstanding:

- how does one send asynchronous javascript URL requests?

- we also need some way to package up images in PyMOL without passaging through the file systems

- we also need to decide what this initial user interface is going to do

"""
