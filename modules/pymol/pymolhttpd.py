# Copyright (C) Schrodinger, LLC.
# All Rights Reserved
#
# For more information, see LICENSE in PyMOL's home directory.
#
# pymolhttpd.py
#
# web server interface for controlling PyMOL

# we make extensive use of Python's build-in in web infrastructure

import sys

if True:
    import http.server as BaseHTTPServer
    import io as StringIO
    import urllib.parse as urlparse
    from urllib.request import urlopen

import cgi
import socket

# we also rely upon Python's json infrastructure

try:
    import simplejson as json
except:
    import json

# standard Python dependencies

import types, os, sys, traceback, threading

# NOTE: Let's attempt to follow Python PEP 8 for coding style for this
# source code file.   URL: http://www.python.org/de/peps/pep-0008
#
# * maximum target line length to be 79 characters.....................seventy9
# * methods and attribute names as lower_case_underscore
# * class names as UpperCaseCaps
# * private symbols start with a leading underscore
# * uniform indentation consisting of 4 spaces (no tabs!)

_json_mime_types = [ 'text/json', 'application/json' ]

class _PymolHTTPRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):

    # for now, we're using a single-threaded server

    # our actual HTTP server class is private for the time being
    # if we need to, then we'll change this

    def wfile_write(self, s):
        if not isinstance(s, bytes):
            s = s.encode('utf-8')
        self.wfile.write(s)

    def do_GET(self):
        self.process_request()

    def do_POST(self):
        self.process_request()

    def log_message(self, format, *args):
        if self.server.pymol_logging:
            BaseHTTPServer.BaseHTTPRequestHandler.log_message(self,format,
                                                              *args)
    def process_request(self):
        """
        parse any URL or FORM arguments and process the request
        """
        # verify that the request is coming from this machine
        try:
            host, port = self.client_address
            if (host[0:6] != '127.0.'):
                self.send_error(403,
                                "Only localhost requests are allowed (not: %s)"
                                % host)
            else:
                self.session = self.server.pymol_session # local session
                self.callback = None
                self.parse_args()
                self.process_urlpath()
        except socket.error:
            traceback.print_exc()
            print("broken pipe")
            pass

    def parse_args(self):
        """
        parses URL arguments into a urlpath (before the ?)
        and a cgiFieldStorage object (args after the ?).
        for example:
        http://localhost:8080/apply/pymol.cmd.color?color=blue&selection=benz
        would yield self.fs.getvalue("color")       as "blue"
        and           self.fs.getvalue("selection") as "benz"
        self.urlpath would be "/apply/pymol.cmd.color"
        """
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

    def process_urlpath(self):
        """
        self.urlpath can be a request for a document, or a
        special request, such as apply or getattr
        """
        parts = self.urlpath.split('/')

        # for example:
        # if http://localhost:8080/apply/pymol.cmd.color?...
        # then parts is ['', 'apply', 'pymol.cmd.color...']
        # or if http://localhost:8080/apply?_json=...
        # then parts is ['', 'apply?_json=...']

        if len(parts) < 2: # then it cannot be a PyMOL request
            self.send_doc() # simple file retrieval
        else: # might be a PyMOL request
            if len(parts) == 2: # no method name or trailing slash -> blank
                parts.append('')
            if (parts[1] == 'apply'): # calling a method
                self.pymol_apply(parts[2])
            elif (parts[1] == 'getattr'): # retrieving a property
                self.pymol_getattr(parts[2])
            elif (parts[1] == 'echo'): # for debugging purposes
                self.send_resp_header(200,'text/plain')
                self.echo_args(parts[2])
            else: # simple file retrieval
                self.send_doc()

    def pymol_getattr(self, attr):
        """
        apply the repr method to the requested attr, but only for
        allowed attributes - those stored in the session dictionary
        """
        key = '/getattr/' + attr;
        if key in self.session:
            try:
                result = repr(self.session[key])
                self.send_json_result(result)
            except:
                self.send_error(500,"Unable to get attribute.")
                self.wfile_write(" %s\n" % attr)
                traceback.print_exc(file=self.wfile)
        else:
            self.send_error(404,"Not a recognized attribute")
            self.wfile_write(" %s is not a recognized attribute\n" % attr)

    def wrap_return(self, result, status="OK", indent=None):
        r = { 'status' : status, 'result' : result }
        if self.server.wrap_natives==1:
            return json.dumps(r, indent=indent)
        else:
            return json.dumps(result, indent=indent)

    def send_json_result(self, result):
        """
        send the mime header and result body.  requests that came from
        XMLHTTPRequest have specified they will accept (expect) json
        formatted results.  other requests will have come from
        ordinary GET or POST requests via links or forms
        """
        if self.callback is not None:
            self.send_resp_header(200,'text/javascript')
            self.wfile_write("%s(%s)"%(self.callback,self.wrap_return(result)))

        else:
            accept_mime = self.headers.get('Accept')
            if accept_mime in _json_mime_types:
                self.send_resp_header(200,accept_mime)
                self.wfile_write(self.wrap_return(result))

            else:
                self.send_resp_header(200,'text/html')
                self.wfile_write("PyMOL's JSON response: <pre>")
                self.wfile_write(self.wrap_return(result,indent=4))
                self.wfile_write("</pre>")

    def send_json_error(self, code, message):
        if self.callback is not None:
            self.send_resp_header(code,'text/javascript')
            self.wfile_write("%s(%s)"%(self.callback,self.wrap_return(message,"ERROR")))
        else:
            accept_mime = self.headers.get('Accept')
            if accept_mime in _json_mime_types:
                self.send_resp_header(code,accept_mime)
                self.wfile_write(self.wrap_return(message,"ERROR"))
            else:
                self.send_resp_header(code,'text/html')
                self.wfile_write("PyMOL's JSON response: <pre>")
                self.wfile_write(self.wrap_return(message,"ERROR",indent=4))
                self.wfile_write("</pre>")

    def send_exception_json(self, code, message):
        fp = StringIO.StringIO()
        traceback.print_exc(file=fp)
        tb = fp.getvalue()
        message = message + tb.split('\n')
        response = json.dumps(message)
        if self.callback is not None:
            self.send_resp_header(code, 'text/javascript')
            self.wfile_write("%s(%s)"%(self.callback,response))
        else:
            accept_mime = self.headers.get('Accept')
            if accept_mime in _json_mime_types:
                self.send_resp_header(code,accept_mime)
                self.wfile_write(response)
            else:
                self.send_resp_header(code,'text/html')
                self.wfile_write("PyMOL's JSON response: <pre>")
                self.wfile_write(json.dumps(json.loads(response),indent=4))
                self.wfile_write("</pre>")

    def pymol_apply(self,method):
        """
        apply the appropriate method held in the session dictionary.
        supply the method arguements in the form of key/value
        """
        args = None
        kwds = None
        query_kwds = {}
        send_multi_result_list = False

        for k in self.fs.keys():
            if k[0:1] == '_': # leading-underscore argument (special handling)
                if k == '_callback':
                    self.callback = self.fs.getfirst(k)
                elif k == '_json': # main path for Javascript API
                    method = json.loads(self.fs.getfirst(k))
                    # [ "my_method", [ arg1, ... ] , { 'key1' : 'val1, ... } ]
                    # or
                    # [ [ "my_met1", [ arg1, ... ], { 'key1' : 'val1, ... } ],
                    #   [ "my_met2", [ arg1, ... ], { 'key1' : 'val1, ... } ] ]
                elif k == '_method': # tentative, not in spec -- may disappear
                    # a method name "my_method"
                    method = json.loads(self.fs.getfirst(k))
                elif k == '_args': # tentative, not in spec -- may disappear
                    args = json.loads(self.fs.getfirst(k))
                elif k == '_kwds': # tentative, not in spec -- may disappear
                    kwds = json.loads(self.fs.getfirst(k))
                # other underscore arguments are ignored (not passed on)
            elif k[0:1] != '_':
                query_kwds[k] = self.fs.getfirst(k)

        blocks = []
        if isinstance(method, str):
            # method is merely a string
            if kwds is None:
                kwds = query_kwds
            if args is None:
                args = ()
            if len(method):
                blocks = [ [ method, args, kwds ] ]
        elif isinstance(method, list) and len(method):
            # method is a list
            if not isinstance(method[0], list):
                blocks = [ method ] # contains just [name, args, kwds]
            else:
                blocks = method
                # contains [ [name, arg, kwds], [name, args, kwds], ... ]
                send_multi_result_list = False # only return final result
        else:
            self.send_json_error(500,[ "Unable to apply method:", str(method)])
            return

        result = []
        if len(blocks):
            for block in blocks:
                if self.server.pymol_logging:
                    print('applying: ' + str(block))
                fn = self.session.get(block[0],None)
                if fn is None and block[0].startswith('pymol.cmd.'):
                    fn = getattr(self.server.pymol_cmd, block[0][10:], None)
                if fn is not None:
                    len_block = len(block)
                    if len_block>1:
                        args = tuple(block[1])
                    else:
                        args = ()
                    if len_block>2:
                        kwds = block[2]
                    else:
                        kwds = {}
                    try:
                        result.append( fn(*args, **kwds) )
                    except:
                        self.send_exception_json(500,
                                                 [ "Exception in: %s" %
                                                   block[0],
                                                   "Args: " + str(args) ,
                                                   "Kwds: " + str(kwds)])
                        return
                else:
                    self.send_json_error(500,[ "Method not found:",
                                               str(block) ])
                    return

                if block[0] == '_quit': # special quit behavior
                    self.send_resp_header()
                    self.wfile_write("<html>")
                    href = None
                    if "href" in kwds:
                        href = str(kwds['href'])
                    elif len(args):
                        href = str(args[1])
                    if href is None:
                        self.wfile_write("<body>")
                    elif not len(href): # simply
                        self.wfile_write("<body onload=\"window.close()\">")
                    else:
                        self.wfile_write(
                            "<body onload=\"document.location.replace('"+
                                         kwds['href']+"')\">")
                    self.wfile_write("<p>PyMOL-HTTPd: Shutting down...</p>")
                    self.wfile_write("<p><i>Please close this window.</i></p>")
                    self.wfile_write("</body></html>")
                    self.wfile.flush()
                    self.server.pymol_cmd.quit()
                    return

        if send_multi_result_list:
            self.send_json_result(result)
        elif len(result):
            self.send_json_result(result[-1])
        else:
            self.send_json_result(None)
        return

    def send_doc(self):
        """
        send a document (file) in the current directory or any sub-directory
        """
        path_list = self.path.split('/')[1:]
        if '..' in path_list: # prevent access to parent directories
            self.send_error(404,"Illegal path.")
            self.wfile_write(": %s" % self.path)
        elif self.server.pymol_root is None:
            self.send_error(404,"No content root specified.")
        else:
            try:
                full_path = os.path.join(*[self.server.pymol_root] +
                                         list(path_list))
                if os.path.isdir(full_path):
                    full_path = full_path + "/index.html"
                fp = open(full_path,"rb")
                self.send_resp_header(200,self.guess_mime(full_path))
                self.wfile_write(fp.read())
                fp.close()
            except:
                self.send_error(404,"Unable to locate document.")
                self.wfile_write(": %s" % self.path)
                self.wfile_write(str(sys.exc_info()))
                # exc_info() is thread safe
                # self.wfile.write(sys.exc_value) # exc_value not thread safe

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

    def send_custom_headers(self):
        for item in getattr(self.server, 'custom_headers', ()):
            self.send_header(item[0], item[1])

    def send_error(self,errcode,errmsg):
        self.send_response(errcode)
        self.send_header('Content-type', 'text/plain')
        self.send_header('Pragma','no-cache')
        self.send_header('Cache-Control','no-cache, must-revalidate')
        self.send_header('Expires','Sat, 10 Jan 2008 01:00:00 GMT')
        self.send_custom_headers()
        self.end_headers()
        self.wfile_write("PyMOL-HTTPd-Error: "+errmsg+"\n")

    def send_resp_header(self, code=200, mime='text/html'):
        self.send_response(code)
        self.send_header('Content-type', mime)
        self.send_header('Pragma','no-cache')
        self.send_header('Cache-Control','no-cache, must-revalidate')
        self.send_header('Expires','Sat, 10 Jan 2008 01:00:00 GMT')
        self.send_custom_headers()
        self.end_headers()

    def echo_args(self):
        """
        for debugging requests
        """
        self.wfile_write("%s\n" % self.command)
        if (self.fs):
            for k in self.fs.keys():
                self.wfile_write("%s = " % k)
                # key can have multiple values, as with checkboxes,
                # but also arbitrarily
                if (isinstance(self.fs[k], list)):
                    self.wfile_write("%s\n" % self.fs.getlist(k))
                else:
                    # key can be uploaded file
                    if (self.fs[k].filename):
                        self.wfile_write("%s\n" % self.fs[k].filename)
                        fp = self.fs[k].file
                        #self.wfile.write("FILE %s" % cgi.escape(repr(fp)))
                        #self.wfile.write("%s\n" % fp.name)
                        # fails for StringIO instances
                        self.wfile_write("%s\n" % repr(fp))
                        # two ways to get file contents
                        #file_contents = self.fs.getvalue(k)
                        #file_contents = fp.read()
                        #self.wfile.write("%s" % file_contents)
                    else:
                        #plain-old key/value
                        self.wfile_write("%s\n" % self.fs.getvalue(k))
        else:
            self.wfile_write("No args\n")

# this is the public class we're exposing to PyMOL consortium members

class PymolHttpd:

    def __init__(self, port=8080, root=None, logging=1, wrap_natives=0, self_cmd=None, headers=()):
        '''
        :param headers: A list or tuple of (key, value) header items to send with each response.
        '''
        if self_cmd is None:
            # fallback on the global singleton PyMOL API
            try:
                from pymol import cmd
                self_cmd = cmd
            except ImportError:
                self_cmd = None
        self.port = int(port)
        self.stop_event = threading.Event()
        self.stop_event.set()
        self.root = root
        self.cmd = self_cmd
        session = {}
        self.session = session

        # Special methods for the web interface

        session['_quit'] = lambda href=None,s=self:s.quit()

        # JavaScript workarounds for keyword clashes

        session['pymol.cmd.delete_'] = self_cmd.delete
        session['pymol.cmd.super_'] = self_cmd.super

        ## Unsafe methods to workaround (uses eval)

        session['pymol.cmd.label'] = self_cmd.label2 # no-eval version

        self.server = BaseHTTPServer.HTTPServer(('', self.port),
                                                _PymolHTTPRequestHandler)
        self.server.wrap_natives = wrap_natives
        self.server.custom_headers = headers

        if self.port == 0:
            self.port = self.server.socket.getsockname()[1]
        self.server.pymol_session = self.session
        self.server.pymol_root = self.root
        if self.root is not None:
            os.environ['PYMOL_HTTP_ROOT'] = self.root
        self.server.pymol_cmd = self.cmd
        self.server.pymol_logging = logging

    def _server_thread(self):
        while not self.stop_event.isSet():
            self.server.handle_request()

    def start(self):
        print ( " PyMOL-HTTPd: serving requests on http://localhost:%d" %
                self.port )
        t = threading.Thread(target=self._server_thread)
        t.setDaemon(1)
        self.stop_event.clear()
        t.start()

    def stop(self):
        if not self.stop_event.isSet():
            self.stop_event.set()
            try: # create a request in order to release the handler
                urlopen("http://localhost:%d" % self.port)
            except:
                pass
            self.server.socket.close()

    def quit(self):
        self.stop_event.set()

    def expose(self, name, value):
        '''
        exposes a Python method or symbol to the web services interface
        '''
        self.session[name] = value

# default behavior if run explicitly from PyMOL

if __name__ == 'pymol': # launched inside PyMOL

    # initialize the server

    server = PymolHttpd()

    # handle_requests (fires off a separate thread)

    server.start()
