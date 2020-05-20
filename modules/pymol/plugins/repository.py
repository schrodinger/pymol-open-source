'''
PyMOL Plugins Engine, Repository Support

This code is in very experimental state!

(c) 2011-2012 Thomas Holder, PyMOL OS Fellow
License: BSD-2-Clause

'''

import sys
if True:
    import urllib.request as urllib2
    from urllib.parse import urlparse
    from urllib.error import URLError, HTTPError

from .installation import supported_extensions

def urlopen(url):
    '''
    urlopen wrapper with timeout from preferences ("network_timeout").

    The timeout does not effect the "urlopen" call itself, that takes 20sec
    for unavailable urls in my tests. Also "socket.setdefaulttimeout" does
    not change that.
    '''
    from . import pref_get

    timeout = pref_get('network_timeout', 10.0)
    return urllib2.urlopen(url, timeout=timeout)


def urlreadstr(url, encoding='iso-8859-1'):
    '''
    Download helper to obtain 'str' content with Python 3.

    If `encoding` is None, then return `bytes`.
    '''
    handle = urlopen(url)
    content = handle.read()

    if encoding:
        charset = handle.headers.get_content_charset() or encoding
        content = content.decode(charset, errors='ignore')

    handle.close()

    return content


class Repository():
    '''
    Abstract repository class

    All open/read operations should raise IOError on failure.
    '''
    def __init__(self, url):
        self.url = url

    def list(self):
        '''
        Return a list of filenames
        '''
        try:
            return self.list_indexfile()
        except:
            return self.list_scan()

    def list_indexfile(self):
        s = self.retrieve('pluginindex.txt', binary=False)
        return s.splitlines()

    def list_scan(self):
        raise NotImplementedError

    def retrieve(self, name, binary=True):
        '''
        Return file content as string
        '''
        raise NotImplementedError

    def copy(self, name, dst):
        '''
        Copy file. The destination may be a directory. Returns the copy file name.
        '''
        import os

        if os.path.isdir(dst):
            dst = os.path.join(dst, os.path.basename(name))

        content = self.retrieve(name)
        f = open(dst, 'wb')
        f.write(content)
        f.close()

        return dst

    def is_supported(self, name):
        if len(name) == 0 or name[0] in ['.', '_']:
            return False
        for ext in supported_extensions:
            if name.endswith(ext):
                return True
        return False

    def get_full_url(self, name):
        import re
        baseurl = re.sub(r'/[^/]*$', '', self.url)
        return baseurl + '/' + name

class HttpRepository(Repository):
    '''
    HTML page over HTTP
    '''
    def list_scan(self):
        import re

        # fetch as string
        content = urlreadstr(self.url)

        # clear comments
        re_comment = re.compile(r'<!\s*--.*?--\s*>', re.DOTALL)
        content = re_comment.sub('', content)

        # find links
        names = []
        re_a = re.compile(r'<a\s+(.*?)>')
        re_href = re.compile(r'''href\s*=\s*(\S+|".+?"|'.+?')''')
        re_anchor = re.compile(r'#.*')
        for attribs in re_a.findall(content):
            for name in re_href.findall(attribs):
                if name[0] in ['"', "'"]:
                    if name[0] != name[-1]:
                        continue
                    name = name[1:-1]
                if '#' in name:
                    name = re_anchor.sub('', name)
                # filter for supported types
                if self.is_supported(name):
                    names.append(name)

        return names

    def retrieve(self, name, binary=True):
        url = self.get_full_url(name)

        if binary:
            handle = urlopen(url)
            content = handle.read()
            handle.close()
        else:
            content = urlreadstr(url)

        return content

    def get_full_url(self, name):
        import re
        if '://' in name:
            return name
        if name.startswith('/'):
            baseurl = '/'.join(self.url.split('/')[:3])
        else:
            baseurl = re.sub(r'/[^/]*$', '', self.url)
        return baseurl + '/' + name

class GithubRepository(HttpRepository):
    '''
    http://developer.github.com/v3/

    Inherit HttpRepository to reuse "retrieve" method.
    '''

    def __init__(self, url):
        import re
        m = re.search(r'github.com[:/]([^/]+)/([^/]+)', url)

        if m is None:
            raise KeyError('cannot parse github.com url')

        self.user = m.group(1)
        self.repo = m.group(2)
        if self.repo.endswith('.git'):
            self.repo = self.repo[:-4]

        # for now: only support main repository
        assert self.user == 'Pymol-Scripts'
        assert self.repo == 'Pymol-script-repo'

        # for HttpRepository.retrieve
        self.url = 'https://github.com/%s/%s/raw/master/' % (self.user, self.repo)

    def list_scan(self):
        sha = 'master'
        r = self.fetchjson('/repos/%s/%s/git/trees/%s' % (self.user, self.repo, sha))

        names = []
        for d in r['tree']:
            if d['type'] != 'blob':
                continue
            name = d['path']
            if self.is_supported(name):
                names.append(name)

        return names

    list = list_scan

    def fetchjson(self, url):
        import json
        content = urlreadstr('https://api.github.com' + url)
        return json.loads(content)

class LocalRepository(Repository):
    def __init__(self, url):
        r = urlparse(url)
        self.url = r.path

    def list_scan(self):
        import os
        names = os.listdir(self.url)
        return list(filter(self.is_supported, names))

    def retrieve(self, name, binary=True):
        url = self.get_full_url(name)
        handle = open(url, "rb" if binary else "r")
        content = handle.read()
        handle.close()
        return content

    def copy(self, name, dst):
        import shutil
        url = self.get_full_url(name)
        shutil.copy(url, dst)
        return dst

    def get_full_url(self, name):
        import os
        return os.path.join(self.url, name)

def guess(url):
    u = url.lower()

    if 'github.com' in u:
        return GithubRepository(url)

    if u.startswith('http://') or \
            u.startswith('https://'):
        return HttpRepository(url)

    if u.startswith('file://') or \
            u.startswith('/'):
        return LocalRepository(url)

    raise KeyError('cannot guess repository type')

def fetchscript(title, dest=None, run=1, quiet=1):
    '''
DESCRIPTION

    Fetch script from PyMOLWiki or Pymol-script-repo

    Returns filename on success, or None on failure.

ARGUMENTS

    title = string: Wiki page title or full URL
    '''
    import re, os
    from pymol import cmd, CmdException

    title = title.strip()
    quiet = int(quiet)
    if dest is None:
        dest = cmd.get('fetch_path')

    # remove hash
    if '#' in title:
        title = title.split('#')[0]

    # github
    git_master = 'https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/'
    m = re.match(r'https://(raw\.)?github\.com/.+', title)
    if m is not None:
        filename = url = title
        if not m.group(1):
            url = url.replace('/blob/', '/raw/', 1)
        rawscript = 1

    # pymolwiki
    elif '://' not in title or 'pymolwiki.org' in title:
        if title.startswith('http'):
            a = title.split('pymolwiki.org/index.php/', 2)
            if len(a) == 2:
                title = a[1]
            else:
                m = re.search(r'[?&]title=([^&]+)', title)
                if m is not None:
                    title = m.group(1)
                else:
                    raise CmdException('Failed to parse URL: ' + title)

        title = title[0].upper() + title[1:].replace(' ','_')
        url = "http://pymolwiki.org/index.php?title=%s&action=raw" % (title)
        filename = title + '.py'
        rawscript = 0

    # any url
    else:
        filename = url = title
        rawscript = 1

    filename = os.path.join(dest, filename.rsplit('/')[-1])

    if os.path.exists(filename):
        if not quiet:
            print('File "%s" exists, will not redownload')
    else:
        if not quiet:
            print('Downloading', url)

        # get page content
        try:
            content = urlreadstr(url, None if rawscript else 'utf-8')
        except IOError as e:
            raise CmdException(e, "Plugin-Error")

        if not rawscript:
            # redirect
            redirect = re.match(r'\s*#REDIRECT\s*\[\[(.*?)\]\]', content)
            if redirect is not None:
                return fetchscript(redirect.group(1), dest, run, quiet)

            # parse Infobox
            pattern2 = re.compile(r'\{\{Infobox script-repo.*?\| *filename *= *([^|\s]+).*?\}\}', re.DOTALL)
            chunks2 = pattern2.findall(content)
            if len(chunks2) > 0:
                try:
                    return fetchscript(git_master + chunks2[0], dest, run, quiet)
                except HTTPError:
                    print('Warning: Infobox filename found, but download failed')

            # parse for <source ...>...</source>
            pattern = re.compile(r'<(?:source|syntaxhighlight)\b[^>]*>(.*?)</(?:source|syntaxhighlight)>', re.DOTALL)
            chunks = pattern.findall(content)

            # check script-chunks for cmd.extend
            chunks = [s for s in chunks if 'cmd.extend' in s]

            if len(chunks) == 0:
                raise CmdException('No <source> or <syntaxhighlight> block with cmd.extend found')
            if len(chunks) > 1:
                print('Warning: %d chunks found, only saving first' % (len(chunks)))

            content = chunks[0]

            with open(filename, 'w') as handle:
                handle.write(content)
        else:
            with open(filename, 'wb') as handle:
                handle.write(content)

    if int(run):
        cmd.do("run " + filename, echo=not quiet)

    return filename

if __name__ == '__main__':
    try:
        import socket_hack
    except ImportError:
        pass

    r1 = guess('http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/')
    print(r1.list())

    r1 = guess('https://github.com/Pymol-Scripts/Pymol-script-repo')
    print(r1.list())

    r1 = guess('/opt/pymol-svn/modules/pmg_tk/startup')
    print(r1.list())

# vi:expandtab:smarttab:sw=4
