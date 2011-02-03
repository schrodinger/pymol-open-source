# Copyright (C) Schrodinger, LLC.
# All Rights Reserved
#
# For more information, see LICENSE in PyMOL's home directory.
#
# justhttpd.py
#
# vanilla web server designed for testing multi-origin applications
# by serving up content on 127.0.0.1:xxxx instead of localhost:yyyy

import BaseHTTPServer, cgi, urlparse, socket

import types, os, sys, traceback, threading

class _HTTPRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):

    def do_GET(self):
        self.process_request()

    def do_POST(self):
        self.process_request()

    def process_request(self):
        """
        parse any URL or FORM arguments and process the request
        """
        # verify that the request is coming from localhost
        try:
            host, port = self.client_address
            if host != '127.0.0.1':
                self.send_error(403,
                                "Only localhost requests are allowed (not: %s)"
                                % host)
            else:
                self.callback = None 
                self.parse_args()
                self.send_doc()
        except socket.error:
            pass
        
    def parse_args(self):
        if (self.command == "POST"):
            self.fs = cgi.FieldStorage(fp=self.rfile, headers=self.headers,
                                       environ = {'REQUEST_METHOD':'POST'},
                                       keep_blank_values = 1)
            self.urlpath = self.path
        elif (self.command == "GET"):
            scheme,netloc,path,params,qs,fragment = urlparse.urlparse(self.path)
            self.fs = cgi.FieldStorage(environ = {'REQUEST_METHOD':'GET',
                                                  'QUERY_STRING':qs},
                                       keep_blank_values = 1)
            self.urlpath = path
        else:
            self.fs = None
        
    def send_doc(self):
        """
        send a document (file) in the current directory or any sub-directory
        """
        path_list = self.path.split('/')[1:]
        if '..' in path_list: # prevent access to parent directories
            self.send_error(404,"Illegal path.")
            self.wfile.write(": %s" % self.path)
        elif self.server.root == None:
            self.send_error(404,"No content root specified.")
        else:
            try:
                full_path = os.path.join(*[self.server.root] +
                                         list(path_list))
                print full_path
                if os.path.isdir(full_path):
                    full_path = full_path + "/index.html"
                fp = open(full_path,"rb")
                self.send_ok(self.guess_mime(full_path))
                self.wfile.write(fp.read())
                fp.close()
            except:
                self.send_error(404,"Unable to locate document.")
                self.wfile.write(": %s" % self.path)
                self.wfile.write(str(sys.exc_info())) # exc_info() is thread safe
                # self.wfile.write(sys.exc_value) # exc_value is not thread safe

    def guess_mime(self,path):
        """
        guess the mime type based on the file extension
        """
        if path.endswith('.html'):
            return 'text/html'
        elif path.endswith('.js'):
            return 'application/x-javascript'
        elif path.endswith('.jpg'):
            return 'image/jpeg'
        elif path.endswith('.png'):
            return 'image/png'
        elif path.endswith('.gif'):
            return 'image/gif'
        elif path.endswith('.sdf'):
            return 'chemical/x-mdl-sdfile'
        elif path.endswith('.mol'):
            return 'chemical/x-mdl-molfile'
        elif path.endswith('.pwg'):
            return 'application/x-pymol'
        else:
            return 'text/plain'
            
    def send_error(self,errcode,errmsg):
        try:
            self.send_response(errcode)
            self.send_header('Content-type', 'text/plain')
            self.end_headers()
            self.wfile.write("HTTPd-Error: "+errmsg+"\n")
        except:
            # right now we're swallowing any/all exceptions
            # (e.g. Broken Pipe)
            pass
        
    def send_ok(self, mime='text/html'):
        self.send_response(200)
        self.send_header('Content-type', mime)
        self.send_header('Pragma','no-cache')
        self.send_header('Cache-Control','no-cache, must-revalidate')
        self.send_header('Expires','Sat, 10 Jan 2008 01:00:00 GMT')
        self.end_headers()

    def echo_args(self):
        """
        for debugging requests
        """
        self.wfile.write("%s\n" % self.command)
        if (self.fs):
            for k in self.fs.keys():
                self.wfile.write("%s = " % k)
                # key can have multiple values, as with checkboxes,
                # but also arbitrarily
                if (isinstance(self.fs[k], types.ListType)):
                    self.wfile.write("%s\n" % self.fs.getlist(k))
                else:
                    # key can be uploaded file
                    if (self.fs[k].filename):
                        self.wfile.write("%s\n" % self.fs[k].filename)
                        fp = self.fs[k].file
                        #self.wfile.write("FILE %s" % cgi.escape(repr(fp)))
                        #self.wfile.write("%s\n" % fp.name)
                        # fails for StringIO instances
                        self.wfile.write("%s\n" % repr(fp))
                        # two ways to get file contents
                        #file_contents = self.fs.getvalue(k)
                        #file_contents = fp.read()
                        #self.wfile.write("%s" % file_contents)
                    else:
                        #plain-old key/value
                        self.wfile.write("%s\n" % self.fs.getvalue(k))
        else:
            self.wfile.write("No args\n")

class PlainHttpd:

    def __init__(self, port=0, root=None):
        self.port = int(port)
        self.stop_event = threading.Event()
        self.stop_event.set()
        self.root = root

        self.server = BaseHTTPServer.HTTPServer(('', self.port),
                                                _HTTPRequestHandler)
        if self.port == 0:
            self.port = self.server.socket.getsockname()[1]

        self.server.root = self.root
        
    def _server_thread(self):
        while not self.stop_event.isSet():
            self.server.handle_request()

    def start(self): # spawn thread
        print " HTTPd: serving requests on http://127.0.0.1:%d" % self.port
        t = threading.Thread(target=self._server_thread)
        t.setDaemon(1)
        self.stop_event.clear()
        t.start()

    def stop(self):
        if not self.stop_event.isSet():
            self.stop_event.set()
            try: # create a request in order to release the handler
                import urllib
                urllib.urlopen("http://localhost:%d" % self.port)
            except:
                pass
            self.server.socket.close()
        
def main():

    import os
    
    # initialize the server, with current local working directory as root
    
    server = PlainHttpd(0, ".")

    # get a dynamically assigned port number

    port = server.port
    
    # start handling requests

    server.start()

    # now launch a browser pointing at our server

    import webbrowser
    webbrowser.open("http://127.0.0.1:%d"%port)

if __name__ in [ '__main__', 'pymol' ]: 

    # intended to be launched with normal Python or
    # pymol -qc justhttpd.py

    main()
    import time
    while 1:
        time.sleep(1)
        
