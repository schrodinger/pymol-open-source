import SocketServer
import BaseHTTPServer

import time
import cgi
import threading
import traceback

from pymol import cmd

_server = None
def _shutdown(self_cmd=cmd):
    if _server != None:
        _server.socket.close()
    self_cmd.quit()

# Note, this handler assumes PyMOL is running as a global singleton

def get_status(out, self_cmd=cmd):
    out.write('<html>\n')
    out.write('<header>\n')
    out.write('<script type="text/javascript" src="pymol.js"></script>\n')
    out.write('</header><body>\n')
    out.write('<h3>PyMOL WebGUI Proof of Concept</h3>\n')
    out.write('<table><tr>\n')
    out.write('<td><form action="./status.pymol"><button type="submit">Refresh</button></form></td>\n')
    out.write('<td><form target="_new" action="./ray.pymol?t=%f"><button type="submit">Ray</button></form></td>\n'%
              time.time())    
    out.write('<td><form target="_new" action="./monitor.pymol?t=%f"><button type="submit">Monitor</button></form></td>\n'%
              time.time())    
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
    out.write('<a href="#" onClick="updateImage()"><img src="./draw.pymol?t=%f">'%time.time()+"</img></a>")
    out.write('</body></html>\n')

def get_monitor(out, self_cmd=cmd):
    out.write('<html>\n')
    out.write('<header>\n')
    out.write('<script type="text/javascript" src="pymol.js"></script>\n')
    out.write('</header><body onload="monitorOnLoad()">\n')
    out.write('<img src="./draw.pymol?t=%f">'%time.time()+"</img>")
    out.write('</body></html>\n')

def get_start(out, self_cmd=cmd):
    window_open="javascript: window.open('./status.pymol','hello', 'location=no,toolbar=no,width=400,height=600,top=0,left=880');"
    out.write('<html><body onload="'+window_open+'">\n')
    out.write('<a href="./start.pymol" onclick="'+window_open+'">Launch WebGUI\n')
    out.write('</body></html>')

def write_image(out, ray=0, self_cmd=cmd):
    if ray:
        self_cmd.ray()
    # encode the file descriptor into the PNG filename
    if self_cmd.png(chr(1)+str(out.fileno()),prior=-1) != 1:
        # no prior image available, so wait for update / finish
        self_cmd.sync()
        
class PymolHandler(BaseHTTPServer.BaseHTTPRequestHandler):


    def log_message(self, format, *args):
        pass
        # nuke logging feature for the time being
    
    def do_js(self):
        self.send_response(200)
        self.send_header('Content-type','text/javascript')
        self.end_headers()
        self.wfile.write('''
            

function updateImage()
{
    images = document.getElementsByTagName("img");

    for( var i = 0; i < images.length; i++ ) {
       images[i].src = "./draw.pymol?t=" + new Date().getTime();
    }
    return false;
}

function monitorOnLoad(event)
{
  setInterval('updateImage()',1000)
}

            ''')
            
    def do_pymol(self):
        if "ray.pymol" in self.path: # send image
            self.send_response(200)
            self.send_header('Content-type',	'image/x-png')
            self.send_header('Cache-control', 'no-cache')
            self.send_header('Pragma', 'no-cache')
            self.end_headers()
            write_image(self.wfile,1)
        elif "draw.pymol" in self.path:
            self.send_response(200)
            self.send_header('Content-type',	'image/x-png')
            self.send_header('Cache-control', 'no-cache')
            self.send_header('Pragma', 'no-cache')
            self.end_headers()
            write_image(self.wfile)
        else:
            if "load" in self.path: # load a structure
                cmd.load("$TUT/1hpv.pdb")
                cmd.rock()
            self.send_response(200)
            self.send_header('Content-type',	'text/html')
            self.end_headers()
            if "status.pymol" in self.path:
                get_status(self.wfile)
            elif "monitor.pymol" in self.path:
                get_monitor(self.wfile)
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
                try:
                    self.do_pymol()
                except:
                    traceback.print_exc()
            elif doc.endswith('.js'): # Javascript
                self.do_js()
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

def main():
    try:
        global _server
        _server = BaseHTTPServer.HTTPServer(('', 8080), PymolHandler)
        print 'started httpserver...'
        _server.serve_forever()
    except KeyboardInterrupt:
        print '^C received, shutting down server'
        _server.socket.close()

def open_browser():
    import webbrowser
    time.sleep(1)
    webbrowser.open('http://localhost:8080/status.pymol')
#    import os
#    os.system('open http://localhost:8080/start.pymol')

if __name__ == '__main__':
    main()

if __name__ == 'pymol':

    cmd.set("image_copy_always") # copy all updates into image buffer

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
