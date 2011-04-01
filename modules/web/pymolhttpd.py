# Copyright (C) Schrodinger, LLC.
# All Rights Reserved
#
# For more information, see LICENSE in PyMOL's home directory.
#
# pymolhttpd.py
#
# web server interface for controlling PyMOL

# we make extensive use of Python's build-in in web infrastructure

import BaseHTTPServer, cgi, urlparse
import StringIO, socket

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
            print "broken pipe"
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
        if self.session.has_key(key):
            try:   
                result = repr(self.session[key])
                self.send_json_result(result)
            except:
                self.send_error(500,"Unable to get attribute.")
                self.wfile.write(" %s\n" % attr)
                traceback.print_exc(file=self.wfile)
        else:
            self.send_error(404,"Not a recognized attribute")
            self.wfile.write(" %s is not a recognized attribute\n" % attr)
        
    def send_json_result(self, result):
        """
        send the mime header and result body.  requests that came from
        XMLHTTPRequest have specified they will accept (expect) json
        formatted results.  other requests will have come from
        ordinary GET or POST requests via links or forms
        """
        if self.callback != None:
            self.send_resp_header(200,'text/javascript')
            self.wfile.write("%s(%s)"%(self.callback,json.dumps(result)))

        else:
            accept_mime = self.headers.getheader('Accept')
            if accept_mime in _json_mime_types:
                self.send_resp_header(200,accept_mime)
                self.wfile.write(json.dumps(result))
            else:
                self.send_resp_header(200,'text/html')
                self.wfile.write("PyMOL's JSON response: <pre>")
                self.wfile.write(json.dumps(result,indent=4))
                self.wfile.write("</pre>")
            
    def send_json_error(self, code, message):
        if self.callback != None:
            self.send_resp_header(code,'text/javascript')
            self.wfile.write("%s(%s)"%(self.callback,json.dumps(message)))
        else:
            accept_mime = self.headers.getheader('Accept')            
            if accept_mime in _json_mime_types:
                self.send_resp_header(code,accept_mime)
                self.wfile.write(json.dumps(message))
            else:
                self.send_resp_header(code,'text/html')
                self.wfile.write("PyMOL's JSON response: <pre>")
                self.wfile.write(json.dumps(message,indent=4))
                self.wfile.write("</pre>")

    def send_exception_json(self, code, message):
        fp = StringIO.StringIO()
        traceback.print_exc(file=fp)
        tb = fp.getvalue()
        message = message + tb.split('\n')
        response = json.dumps(message)
        if self.callback != None:
            self.send_resp_header(code, 'text/javascript')
            self.wfile.write("%s(%s)"%(self.callback,response))
        else:
            accept_mime = self.headers.getheader('Accept')
            if accept_mime in _json_mime_types:
                self.send_resp_header(code,accept_mime)
                self.wfile.write(response)
            else:
                self.send_resp_header(code,'text/html')
                self.wfile.write("PyMOL's JSON response: <pre>")
                self.wfile.write(json.dumps(json.loads(response),indent=4))
                self.wfile.write("</pre>")

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
        if isinstance(method,types.StringType):
            # method is merely a string 
            if kwds == None:
                kwds = query_kwds
            if args == None:
                args = ()
            if len(method):
                blocks = [ [ method, args, kwds ] ]
        elif isinstance(method,types.ListType) and len(method):
            # method is a list
            if not isinstance(method[0],types.ListType):
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
                    print 'applying: ' + str(block)
                fn = self.session.get(block[0],None)
                if fn != None:
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
                    self.wfile.write("<html>")
                    href = None
                    if kwds.has_key("href"):
                        href = str(kwds['href'])
                    elif len(args):
                        href = str(args[1])
                    if href == None:
                        self.wfile.write("<body>")
                    elif not len(href): # simply 
                        self.wfile.write("<body onload=\"window.close()\">")
                    else:
                        self.wfile.write(
                            "<body onload=\"document.location.replace('"+
                                         kwds['href']+"')\">")
                    self.wfile.write("<p>PyMOL-HTTPd: Shutting down...</p>")
                    self.wfile.write("<p><i>Please close this window.</i></p>")
                    self.wfile.write("</body></html>")
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
            self.wfile.write(": %s" % self.path)
        elif self.server.pymol_root == None:
            self.send_error(404,"No content root specified.")
        else:
            try:
                full_path = os.path.join(*[self.server.pymol_root] +
                                         list(path_list))
                if os.path.isdir(full_path):
                    full_path = full_path + "/index.html"
                fp = open(full_path,"rb")
                self.send_resp_header(200,self.guess_mime(full_path))
                self.wfile.write(fp.read())
                fp.close()
            except:
                self.send_error(404,"Unable to locate document.")
                self.wfile.write(": %s" % self.path)
                self.wfile.write(str(sys.exc_info()))
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
            
    def send_error(self,errcode,errmsg):
        self.send_response(errcode)
        self.send_header('Content-type', 'text/plain')
        self.send_header('Pragma','no-cache')
        self.send_header('Cache-Control','no-cache, must-revalidate')
        self.send_header('Expires','Sat, 10 Jan 2008 01:00:00 GMT')
        self.end_headers()
        self.wfile.write("PyMOL-HTTPd-Error: "+errmsg+"\n")
        
    def send_resp_header(self, code=200, mime='text/html'):
        self.send_response(code)
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

# this is the public class we're exposing to PyMOL consortium members

class PymolHttpd:

    def __init__(self, port=8080, root=None, logging=1, self_cmd=None):
        if self_cmd == None:
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

        ## session['pymol.cmd.label'] = self_cmd.label

        session['pymol.cmd.label'] = self_cmd.label2 # no-eval version
        
        # BEGIN MACHINE-GENERATED CODE

        session['pymol.cmd.load'] = self_cmd.load
        session['pymol.cmd.load_traj'] = self_cmd.load_traj
        session['pymol.cmd.load_png'] = self_cmd.load_png
        session['pymol.cmd.fragment'] = self_cmd.fragment
        session['pymol.cmd.fetch'] = self_cmd.fetch
        session['pymol.cmd.read_mmodstr'] = self_cmd.read_mmodstr
        session['pymol.cmd.read_molstr'] = self_cmd.read_molstr
        session['pymol.cmd.read_sdfstr'] = self_cmd.read_sdfstr
        session['pymol.cmd.read_pdbstr'] = self_cmd.read_pdbstr
        session['pymol.cmd.read_xplorstr'] = self_cmd.read_xplorstr
        session['pymol.cmd.get_pdbstr'] = self_cmd.get_pdbstr
        session['pymol.cmd.get_fastastr'] = self_cmd.get_fastastr
        session['pymol.cmd.copy'] = self_cmd.copy
        session['pymol.cmd.create'] = self_cmd.create
        session['pymol.cmd.extract'] = self_cmd.extract
        session['pymol.cmd.split_states'] = self_cmd.split_states
        session['pymol.cmd.symexp'] = self_cmd.symexp
        session['pymol.cmd.ramp_new'] = self_cmd.ramp_new
        session['pymol.cmd.set_name'] = self_cmd.set_name
        session['pymol.cmd.map_new'] = self_cmd.map_new
        session['pymol.cmd.map_set'] = self_cmd.map_set
        session['pymol.cmd.map_set_border'] = self_cmd.map_set_border
        session['pymol.cmd.map_double'] = self_cmd.map_double
        session['pymol.cmd.map_halve'] = self_cmd.map_halve
        session['pymol.cmd.map_trim'] = self_cmd.map_trim
        session['pymol.cmd.isodot'] = self_cmd.isodot
        session['pymol.cmd.isolevel'] = self_cmd.isolevel
        session['pymol.cmd.isomesh'] = self_cmd.isomesh
        session['pymol.cmd.isosurface'] = self_cmd.isosurface
        session['pymol.cmd.slice_new'] = self_cmd.slice_new
        session['pymol.cmd.gradient'] = self_cmd.gradient
        session['pymol.cmd.ungroup'] = self_cmd.ungroup
        session['pymol.cmd.group'] = self_cmd.group
        session['pymol.cmd.pseudoatom'] = self_cmd.pseudoatom
        session['pymol.cmd.fab'] = self_cmd.fab
        session['pymol.cmd.enable'] = self_cmd.enable
        session['pymol.cmd.disable'] = self_cmd.disable
        session['pymol.cmd.delete'] = self_cmd.delete
        session['pymol.cmd.reinitialize'] = self_cmd.reinitialize
        session['pymol.cmd.deselect'] = self_cmd.deselect
        session['pymol.cmd.select'] = self_cmd.select
        session['pymol.cmd.indicate'] = self_cmd.indicate
        session['pymol.cmd.select_list'] = self_cmd.select_list
        session['pymol.cmd.pop'] = self_cmd.pop
        session['pymol.cmd.angle'] = self_cmd.angle
        session['pymol.cmd.dihedral'] = self_cmd.dihedral
        session['pymol.cmd.dist'] = self_cmd.dist
        session['pymol.cmd.distance'] = self_cmd.distance
        session['pymol.cmd.get_angle'] = self_cmd.get_angle
        session['pymol.cmd.get_dihedral'] = self_cmd.get_dihedral
        session['pymol.cmd.get_distance'] = self_cmd.get_distance
        session['pymol.cmd.get_area'] = self_cmd.get_area
        session['pymol.cmd.color'] = self_cmd.color
        session['pymol.cmd.bg_color'] = self_cmd.bg_color
        session['pymol.cmd.rebuild'] = self_cmd.rebuild
        session['pymol.cmd.refresh'] = self_cmd.refresh
        session['pymol.cmd.recolor'] = self_cmd.recolor
        session['pymol.cmd.set_color'] = self_cmd.set_color
        session['pymol.cmd.set_object_color'] = self_cmd.set_object_color
        session['pymol.cmd.show'] = self_cmd.show
        session['pymol.cmd.show_as'] = self_cmd.show_as
        session['pymol.cmd.hide'] = self_cmd.hide
        session['pymol.cmd.cartoon'] = self_cmd.cartoon
        session['pymol.cmd.spectrum'] = self_cmd.spectrum
        session['pymol.cmd.center'] = self_cmd.center
        session['pymol.cmd.zoom'] = self_cmd.zoom
        session['pymol.cmd.reset'] = self_cmd.reset
        session['pymol.cmd.clip'] = self_cmd.clip
        session['pymol.cmd.orient'] = self_cmd.orient
        session['pymol.cmd.origin'] = self_cmd.origin
        session['pymol.cmd.set_view'] = self_cmd.set_view
        session['pymol.cmd.get_view'] = self_cmd.get_view
        session['pymol.cmd.move'] = self_cmd.move
        session['pymol.cmd.turn'] = self_cmd.turn
        session['pymol.cmd.rock'] = self_cmd.rock
        session['pymol.cmd.stereo'] = self_cmd.stereo
        session['pymol.cmd.get'] = self_cmd.get
        session['pymol.cmd.set'] = self_cmd.set
        session['pymol.cmd.set_bond'] = self_cmd.set_bond
        session['pymol.cmd.unset'] = self_cmd.unset
        session['pymol.cmd.unset_bond'] = self_cmd.unset_bond
        session['pymol.cmd.get_setting_boolean'] = self_cmd.get_setting_boolean
        session['pymol.cmd.get_setting_int'] = self_cmd.get_setting_int
        session['pymol.cmd.get_setting_float'] = self_cmd.get_setting_float
        session['pymol.cmd.get_setting_legacy'] = self_cmd.get_setting_legacy
        session['pymol.cmd.get_setting_tuple'] = self_cmd.get_setting_tuple
        session['pymol.cmd.get_setting_text'] = self_cmd.get_setting_text
        session['pymol.cmd.window'] = self_cmd.window
        session['pymol.cmd.viewport'] = self_cmd.viewport
        session['pymol.cmd.full_screen'] = self_cmd.full_screen
        session['pymol.cmd.quit'] = self_cmd.quit
        session['pymol.cmd.draw'] = self_cmd.draw
        session['pymol.cmd.ray'] = self_cmd.ray
        session['pymol.cmd.align'] = self_cmd.align
        session['pymol.cmd.super'] = self_cmd.super
        session['pymol.cmd.fit'] = self_cmd.fit
        session['pymol.cmd.rms'] = self_cmd.rms
        session['pymol.cmd.rms_cur'] = self_cmd.rms_cur
        session['pymol.cmd.intra_fit'] = self_cmd.intra_fit
        session['pymol.cmd.intra_rms'] = self_cmd.intra_rms
        session['pymol.cmd.intra_rms_cur'] = self_cmd.intra_rms_cur
        session['pymol.cmd.pair_fit'] = self_cmd.pair_fit
        session['pymol.cmd.space'] = self_cmd.space
        session['pymol.cmd.order'] = self_cmd.order
        session['pymol.cmd.edit_mode'] = self_cmd.edit_mode
        session['pymol.cmd.button'] = self_cmd.button
        session['pymol.cmd.config_mouse'] = self_cmd.config_mouse
        session['pymol.cmd.mouse'] = self_cmd.mouse
        session['pymol.cmd.mask'] = self_cmd.mask
        session['pymol.cmd.unmask'] = self_cmd.unmask
        session['pymol.cmd.count_atoms'] = self_cmd.count_atoms
        session['pymol.cmd.get_chains'] = self_cmd.get_chains
        session['pymol.cmd.get_color_index'] = self_cmd.get_color_index
        session['pymol.cmd.get_color_indices'] = self_cmd.get_color_indices
        session['pymol.cmd.get_object_color_index'] = self_cmd.get_object_color_index
        session['pymol.cmd.get_object_list'] = self_cmd.get_object_list
        session['pymol.cmd.get_color_tuple'] = self_cmd.get_color_tuple
        session['pymol.cmd.get_atom_coords'] = self_cmd.get_atom_coords
        session['pymol.cmd.get_extent'] = self_cmd.get_extent
        session['pymol.cmd.get_names'] = self_cmd.get_names
        session['pymol.cmd.get_names_of_type'] = self_cmd.get_names_of_type
        session['pymol.cmd.get_legal_name'] = self_cmd.get_legal_name
        session['pymol.cmd.get_unused_name'] = self_cmd.get_unused_name
        session['pymol.cmd.get_object_matrix'] = self_cmd.get_object_matrix
        session['pymol.cmd.get_phipsi'] = self_cmd.get_phipsi
        session['pymol.cmd.get_position'] = self_cmd.get_position
        session['pymol.cmd.get_raw_alignment'] = self_cmd.get_raw_alignment
        session['pymol.cmd.get_renderer'] = self_cmd.get_renderer
        session['pymol.cmd.get_symmetry'] = self_cmd.get_symmetry
        session['pymol.cmd.get_title'] = self_cmd.get_title
        session['pymol.cmd.get_type'] = self_cmd.get_type
        session['pymol.cmd.get_version'] = self_cmd.get_version
        session['pymol.cmd.id_atom'] = self_cmd.id_atom
        session['pymol.cmd.identify'] = self_cmd.identify
        session['pymol.cmd.index'] = self_cmd.index
        session['pymol.cmd.phi_psi'] = self_cmd.phi_psi
        session['pymol.cmd.matrix_copy'] = self_cmd.matrix_copy
        session['pymol.cmd.matrix_reset'] = self_cmd.matrix_reset
        session['pymol.cmd.rotate'] = self_cmd.rotate
        session['pymol.cmd.translate'] = self_cmd.translate
        session['pymol.cmd.set_object_ttt'] = self_cmd.set_object_ttt
        session['pymol.cmd.set_dihedral'] = self_cmd.set_dihedral
        session['pymol.cmd.transform_object'] = self_cmd.transform_object
        session['pymol.cmd.transform_selection'] = self_cmd.transform_selection
        session['pymol.cmd.translate_atom'] = self_cmd.translate_atom
        session['pymol.cmd.update'] = self_cmd.update
        session['pymol.cmd.attach'] = self_cmd.attach
        session['pymol.cmd.bond'] = self_cmd.bond
        session['pymol.cmd.unbond'] = self_cmd.unbond
        session['pymol.cmd.cycle_valence'] = self_cmd.cycle_valence
        session['pymol.cmd.drag'] = self_cmd.drag
        session['pymol.cmd.dss'] = self_cmd.dss
        session['pymol.cmd.edit'] = self_cmd.edit
        session['pymol.cmd.unpick'] = self_cmd.unpick
        session['pymol.cmd.fix_chemistry'] = self_cmd.fix_chemistry
        session['pymol.cmd.flag'] = self_cmd.flag
        session['pymol.cmd.fuse'] = self_cmd.fuse
        session['pymol.cmd.get_editor_scheme'] = self_cmd.get_editor_scheme
        session['pymol.cmd.h_add'] = self_cmd.h_add
        session['pymol.cmd.h_fill'] = self_cmd.h_fill
        session['pymol.cmd.h_fix'] = self_cmd.h_fix
        session['pymol.cmd.invert'] = self_cmd.invert
        session['pymol.cmd.torsion'] = self_cmd.torsion
        session['pymol.cmd.valence'] = self_cmd.valence
        session['pymol.cmd.clean'] = self_cmd.clean
        session['pymol.cmd.deprotect'] = self_cmd.deprotect
        session['pymol.cmd.protect'] = self_cmd.protect
        session['pymol.cmd.reference'] = self_cmd.reference
        session['pymol.cmd.remove'] = self_cmd.remove
        session['pymol.cmd.remove_picked'] = self_cmd.remove_picked
        session['pymol.cmd.rename'] = self_cmd.rename
        session['pymol.cmd.replace'] = self_cmd.replace
        session['pymol.cmd.sculpt_purge'] = self_cmd.sculpt_purge
        session['pymol.cmd.sculpt_deactivate'] = self_cmd.sculpt_deactivate
        session['pymol.cmd.sculpt_activate'] = self_cmd.sculpt_activate
        session['pymol.cmd.sculpt_iterate'] = self_cmd.sculpt_iterate
        session['pymol.cmd.set_geometry'] = self_cmd.set_geometry
        session['pymol.cmd.set_symmetry'] = self_cmd.set_symmetry
        session['pymol.cmd.set_title'] = self_cmd.set_title
        session['pymol.cmd.smooth'] = self_cmd.smooth
        session['pymol.cmd.sort'] = self_cmd.sort
        session['pymol.cmd.undo'] = self_cmd.undo
        session['pymol.cmd.push_undo'] = self_cmd.push_undo
        session['pymol.cmd.redo'] = self_cmd.redo
        session['pymol.cmd.wizard'] = self_cmd.wizard
        session['pymol.cmd.replace_wizard'] = self_cmd.replace_wizard
        session['pymol.cmd.count_frames'] = self_cmd.count_frames
        session['pymol.cmd.count_states'] = self_cmd.count_states
        session['pymol.cmd.mset'] = self_cmd.mset
        session['pymol.cmd.madd'] = self_cmd.madd
        session['pymol.cmd.mclear'] = self_cmd.mclear
        session['pymol.cmd.mmatrix'] = self_cmd.mmatrix
        session['pymol.cmd.mdump'] = self_cmd.mdump
        session['pymol.cmd.mview'] = self_cmd.mview
        session['pymol.cmd.forward'] = self_cmd.forward
        session['pymol.cmd.backward'] = self_cmd.backward
        session['pymol.cmd.rewind'] = self_cmd.rewind
        session['pymol.cmd.middle'] = self_cmd.middle
        session['pymol.cmd.ending'] = self_cmd.ending
        session['pymol.cmd.mplay'] = self_cmd.mplay
        session['pymol.cmd.mtoggle'] = self_cmd.mtoggle
        session['pymol.cmd.mstop'] = self_cmd.mstop
        session['pymol.cmd.frame'] = self_cmd.frame
        session['pymol.cmd.get_movie_playing'] = self_cmd.get_movie_playing
        session['pymol.cmd.get_state'] = self_cmd.get_state
        session['pymol.cmd.get_frame'] = self_cmd.get_frame
        session['pymol.cmd.view'] = self_cmd.view
        session['pymol.cmd.get_scene_dict'] = self_cmd.get_scene_dict
        session['pymol.cmd.get_scene_list'] = self_cmd.get_scene_list
        session['pymol.cmd.scene'] = self_cmd.scene
        session['pymol.cmd.scene_order'] = self_cmd.scene_order

        # END MACHINE-GENERATED CODE

        self.server = BaseHTTPServer.HTTPServer(('', self.port),
                                                _PymolHTTPRequestHandler)
        if self.port == 0:
            self.port = self.server.socket.getsockname()[1]
        self.server.pymol_session = self.session
        self.server.pymol_root = self.root
        if self.root != None:
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
                import urllib
                urllib.urlopen("http://localhost:%d" % self.port)
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

  
