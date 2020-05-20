import re
import sys

# unicode can't be memoryviewed, so just slice
buffer = lambda s, i=0: s[i:]

from pymol import cmd, parsing

DEBUG = False

STRICT = [parsing.STRICT, parsing.ABORT, parsing.SECURE, parsing.LEGACY]
LITERAL = [
    parsing.LITERAL, parsing.LITERAL1, parsing.LITERAL2, parsing.MOVIE,
    parsing.PYTHON
]
MULTILINE = [parsing.PYTHON_BLOCK, parsing.EMBED, parsing.SKIP]

brackopen = {'(': ')', '[': ']', '{': '}'}


def fancysplit(s, delim=',', stop=(';', '\n'), maxsplit=-1):
    stack = []
    quote = ''
    tristart = -1
    escaped = False

    for i, c in enumerate(s):
        if i == maxsplit:
            break

        if escaped:
            escaped = False
        elif c == '\\':
            escaped = True
        elif quote:
            if c == quote:
                if tristart == -1:
                    quote = ''
                elif i > tristart + 2 and s[i + 1:i + 3] == c + c:
                    tristart = -1
                    quote = ''
        elif c in ('"', "'"):
            quote = c
            if s[i + 1:i + 3] == c + c:
                tristart = i
        elif c in brackopen:
            stack.append(brackopen[c])
        elif stack:
            if c == stack[-1]:
                stack.pop()
        elif c == delim:
            yield i
        elif c in stop:
            yield i
            return

    if stack and DEBUG:
        print('warning: unmatched quotes')

    yield len(s)


def findnewline(s):
    m = re.search(r'(?:^|[^\\])$', s, re.M)
    if m is None:
        return len(s)
    return m.end() + 1


def parse_pml(text):
    i = 0
    L = len(text)

    while i < L:
        expr_start = i

        m = re.match(r'(_ )?[ \t]*($|#|/|@|[^\s;]+)', buffer(text, i), re.M)
        if m is None:
            break

        com = m.group(2)

        hit = {
            'type': 'python',
            'quiet': m.group(1) is not None,
            'command': (m.start(2) + i, m.end(2) + i),
            'multiline': 0,  # 0, 1 (closed), 2 (open)
        }

        mode = parsing.PYTHON

        if com == '':
            i += m.end() + 1
            continue
        elif com == '/':
            pass
        elif com == '#':
            # TODO could be handled as command
            hit['type'] = 'comment'
        elif com == '@':
            # TODO could be handled as command
            hit['type'] = 'command'
            mode = parsing.STRICT
        else:
            com = cmd.kwhash.shortcut.get(com, com)
            if com in cmd.keyword:
                mode = cmd.keyword[com][4]

            if mode != parsing.PYTHON:
                hit['type'] = 'command'
            else:
                i += m.start() - m.end()

        i += m.end()
        arg_i = [i]

        if mode in LITERAL:
            i += findnewline(buffer(text, i))
            arg_i.append(i)
        else:
            arg_i.extend((i + p) for p in fancysplit(buffer(text, i)))
            i = arg_i[-1] + 1

            if mode in MULTILINE:
                m = re.search(r'^(_ )?[ \t]*' + com + r' end\b',
                              buffer(text, arg_i[-1]), re.M)
                if m is None:
                    hit['multiline'] = 2
                    i = len(text)
                    arg_i.append(i)
                else:
                    hit['multiline'] = 1
                    arg_i.append(i + m.start())
                    i += m.end()

        hit['args'] = arg_i
        hit['expr'] = expr_start, i

        yield hit
