# Another simple example of how PyMOL can be controlled using a web browser
# using Python's built-in web server capabilities

try:
    import BaseHTTPServer
except ImportError:
    import http.server as BaseHTTPServer

import time
import cgi
import threading
import traceback
import os, sys, re

from pymol import cmd
from chempy.sdf import SDF

# example 3D sd file

input_sdf = os.environ['PYMOL_PATH']+"/test/dat/ligs3d.sdf"

class SafeDict(dict):
    '''
    we will need to synchronize access if we later adopt a
    multi-threaded approach
    '''
    pass

# global singleton class for holding server state

ServerState = SafeDict()

# we're not maintaining session tokens yet...

default_token = 0

def write_table(out, table):
    out.write('<html>\n')
    out.write('<header>\n')
    out.write('<link rel="stylesheet" type="text/css" href="/pymol.css"></link>')
    out.write('<script type="text/javascript" src="pymol.js"></script>\n')
    out.write('</header>')
    out.write('<body>\n')
    out.write('<form action="./quit.pymol"><button type="submit">Quit</button></form>\n')     
    out.write('<table>\n')
    header = table['header']
    for heading in header:
        out.write('<th>')
        out.write(heading)
        out.write('</th>')
    
    body = table['body']
    for row in body:
        out.write('<tr>')
        for col in row:
            out.write('<td>')
            out.write(col)
            out.write('</td>')
        out.write('</tr>\n')
    out.write('</table>')
    out.write('<iframe name="myStatus" width="500" height="60" frameborder=0>')
    out.write('</iframe>\n')
    out.write('<form name="hidden" target="myStatus"></form>')
    out.write('</body></html>\n')

_server = None
def _shutdown(self_cmd=cmd):
    global _server
    if _server != None:
        _server.socket.close()
    self_cmd.quit()

class PymolHandler(BaseHTTPServer.BaseHTTPRequestHandler):

    def log_message(self, format, *args):
        # nuke logging for the time being since it slows down PyMOL
        pass

    def do_css(self):
        self.send_response(200)
        self.send_header('Content-type','text/css')
        self.end_headers()
        self.wfile.write('''
table {
	border-width: 1px 1px 1px 1px;
	border-spacing: 1px;
	border-style: outset outset outset outset;
	border-color: gray gray gray gray;
	border-collapse: separate;
	background-color: gray;
}
table th {
	border-width: 1px 1px 1px 1px;
	padding: 2px 5px 2px 5px;
	border-style: inset inset inset inset;
	border-color: gray gray gray gray;
	background-color: lightblue;
	-moz-border-radius: 0px 0px 0px 0px;
}
table td {
	border-width: 1px 1px 1px 1px;
	padding: 2px 5px 2px 5px;
	border-style: inset inset inset inset;
	border-color: gray gray gray gray;
	background-color: white;
	-moz-border-radius: 0px 0px 0px 0px;
}

a:link {text-decoration: none; color: blue; }
a:visited {text-decoration: none; color: blue; }
a:active {text-decoration: none; color: blue; }

        ''')
        
    def do_js(self):
        self.send_response(200)
        self.send_header('Content-type','text/javascript')
        self.end_headers()

        self.wfile.write('''

function load(molid)
{
// unnecessary...but keeping it around for later...
// this doesnt actually work anyway!
                         
    document.forms['hidden'].method='get';
    document.forms['hidden'].action='http://localhost:8080/load.pymol?test=1&molid=' + molid;
    document.forms['hidden'].submit();
}

            ''')

    def do_pymol(self):
        if "table.pymol" in self.path: # send table
            session = ServerState[default_token]
            self.send_response(200)
            self.send_header('Content-type',	'text/html')
            self.send_header('Cache-control', 'no-cache')
            self.send_header('Pragma', 'no-cache')
            self.end_headers()
            if 'table' not in session:
                self.wfile.write("<p>No table defined.</p>\n")
            else:
                write_table(self.wfile, session['table'])
        elif "quit.pymol" in self.path:
            self.wfile.write('<html><body><p>Quitting...</p></body></html>')
            self.wfile.flush()
            _shutdown()
        elif "load.pymol" in self.path:
            self.send_response(200)
            self.send_header('Content-type',	'text/html')
            self.send_header('Cache-control', 'no-cache')
            self.send_header('Pragma', 'no-cache')
            self.end_headers()
            mo = re.search("molid\=([A-Z0-9]+)",self.path)
            if mo:
                mol_id = mo.groups(1)[0]
                session = ServerState[default_token]
                mol_dict = session['data']['mol_dict']
                self_cmd = session['cmd']
                if mol_id in self_cmd.get_names('objects'):
                    if mol_id in self_cmd.get_names('objects',enabled_only=1):
                        self.wfile.write("<p>Disabling %s...</p>"%mol_id)
                        self_cmd.disable(mol_id)
                    else:
                        self.wfile.write("<p>Enabling %s...</p>"%mol_id)
                        self_cmd.enable(mol_id)
                else:
                    self.wfile.write("<p>Loading %s...</p>"%mol_id)
                    self_cmd.read_molstr(mol_dict[mol_id], mol_id)
                    self_cmd.show_as("sticks",mol_id)
            else:
                self.wfile.write("<p>Error processing query: %s</p>"%self.path)
        else: # start page
            self.send_response(200)
            self.send_header('Content-type',	'text/html')
            self.send_header('Cache-control', 'no-cache')
            self.send_header('Pragma', 'no-cache')
            self.end_headers()
            self.wfile.write("<p>Unhandled PyMOL request</p>")
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
            elif doc.endswith('.css'): # Javascript
                self.do_css()
            elif doc.endswith('.html'):
                f = open('.'+self.path) # UNSAFE!!!
                self.send_response(200)
                self.send_header('Content-type',	'text/html')
                self.end_headers()
                self.wfile.write(f.read())
                f.close()
        except IOError:
            self.send_error(404,'File Not Found: %s' % self.path)

    def do_POST(self): # not currently used
        global rootnode
        try:
            ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
            if ctype == 'multipart/form-data':
                query=cgi.parse_multipart(self.rfile, pdict)
            self.send_response(301)
            
            self.end_headers()
            upfilecontent = query.get('upfile')
            print("filecontent", upfilecontent[0])
            self.wfile.write('<HTML>POST OK.<BR><BR>');
            self.wfile.write(upfilecontent[0]);
            
        except :
            pass

def table_from_data(data):

    # pull MOLID to the far left
    
    col_id_list = ['MOLID'] + [x for x in data['col_id_list'] if x!='MOLID']
    content = data['content']

    # create the header fields
    
    header = []
    
    for col_id in  col_id_list:
        header.append(col_id)

    # create the body
    
    body = []
    for row_id in data['row_id_list']:
        row = []
        for col_id in col_id_list:
            if col_id == 'MOLID':
                text = content.get( (row_id,col_id),'')
                row.append('<a target="myStatus" href="load.pymol?molid=%s">'%text +
                           text + '</a>')
            else:
                row.append( content.get( (row_id,col_id),'' ))
        body.append(row)
        
    return {
        'header' : header,
        'body' : body,
        }

def data_from_sdf(sdf_file_path):

    mol_dict = {}

    row_id_list = []
    row_id_dict = {}
    
    col_id_list = []
    col_id_dict = {}

    # first pass, load the identifiers, MOL files, and tag names

    col = 0
    sdf = SDF(sdf_file_path)
    while 1:
        rec = sdf.read()
        if not rec:
            break
        mol = rec.get('MOL')

        # get the unique identifier
        mol_id = mol[0].strip()

        # store the MOL record
        mol_dict[mol_id] = ''.join(mol)
        
        # add row (assuming mol_id is unique)
        row_id_list.append(mol_id)
        row_id_dict[mol_id] = None
        
        # add column (if new)
        for key in rec.kees:
            if key != 'MOL':
                if key not in col_id_dict:
                    col_id_list.append(key)
                    col_id_dict[key] = None
        
    # second pass, read the actual data into the table structure

    content = {}
    
    sdf = SDF(sdf_file_path)
    while 1:
        rec = sdf.read()
        if not rec:
            break
        mol_id = rec.get('MOL')[0].strip()

        for key in rec.kees:
            if key != 'MOL':
                content[ (mol_id,key) ] = rec.get(key)[0].strip()

    return {
        'content' : content,
        'row_id_list' : row_id_list,
        'col_id_list' : col_id_list,
        'mol_dict' : mol_dict
        }

def open_browser():
    import webbrowser
    time.sleep(1)
    webbrowser.open('http://localhost:8080/table.pymol')
#    import os
#    os.system('open http://localhost:8080/start.pymol')

def main():
    try:
        global _server
        _server = BaseHTTPServer.HTTPServer(('', 8080), PymolHandler)
        print('started httpserver...')
        _server.serve_forever()
    except KeyboardInterrupt:
        print('^C received, shutting down server')
        _server.socket.close()

if __name__ == '__main__':
    print("this script must be run from within PyMOL")
    
if __name__ == 'pymol':

    
    session = SafeDict()
    session['cmd'] = cmd #  could replaced by instance instead of module
    session['data'] = data_from_sdf(input_sdf)
    session['table'] = table_from_data(session['data'])
    
    ServerState[default_token] = session
    
    
    t = threading.Thread(target=main)
    t.setDaemon(1)
    t.start()

    t = threading.Thread(target=open_browser)
    t.setDaemon(1)
    t.start()

    
"""


"""
