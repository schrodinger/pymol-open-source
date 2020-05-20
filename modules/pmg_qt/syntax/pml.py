from . import *

from .python import highlightPython

STYLES = {
    'comment': textformat('blue'),
    'keyword': textformat('brown', 'bold'),
    'argument': textformat('green'),
    'error': textformat('white', 'bg:red'),
    'skipped': textformat('gray'),
}

MULTILINE_PYTHON = 0x01000
MULTILINE_EMBED  = 0x02000
MULTILINE_SKIP   = 0x04000
MULTILINE        = 0x0F000

CONTINUED_PYTHON = 0x10000
CONTINUED        = 0xF0000

multiline_command = {
    MULTILINE_PYTHON: 'python',
    MULTILINE_EMBED: 'embed',
    MULTILINE_SKIP: 'skip',
}

multiline_state = dict((v, k) for (k, v) in multiline_command.items())


class Highlighter(QtGui.QSyntaxHighlighter):
    """
    Syntax highlighter for the PyMOL command language.
    """

    def highlightBlock(self, text):
        import re
        from pymol import cmd
        from .pmlparser import parse_pml

        state = max(0, self.previousBlockState())
        state_python = -1

        if state & MULTILINE:
            com = multiline_command[state & MULTILINE]
            if re.match(r'(_ )?[ \t]*' + com + r' end\b', text) is not None:
                self.setFormat(0, len(text), STYLES['keyword'])
                self.setCurrentBlockState(0)
            elif state & MULTILINE_PYTHON:
                state = highlightPython(self, state & ~MULTILINE, text)
                self.setCurrentBlockState(state | MULTILINE_PYTHON)
            else:
                self.setFormat(0, len(text), STYLES['skipped'])
                self.setCurrentBlockState(state)

            return

        if state & CONTINUED:
            if state & CONTINUED_PYTHON:
                state = highlightPython(self, state & ~CONTINUED, text)
                state |= CONTINUED_PYTHON

            if not text.endswith('\\'):
                state = 0

            self.setCurrentBlockState(state)
            return

        for hit in parse_pml(text):
            self.setCurrentBlockState(0)
            state_python = -1

            if hit['type'] == 'comment':
                i, j = hit['expr']
                self.setFormat(i, j - i, STYLES['comment'])
            elif hit['type'] == 'python':
                i = hit['args'][0]
                j = hit['expr'][1]
                if j > i and text[min(j, len(text)) - 1] == ':':
                    # blocks not supported
                    self.setFormat(i, j - i, STYLES['error'])
                else:
                    state_python = highlightPython(self, 0, text, i, j)
            elif hit['type'] == 'command':
                i, j = hit['command']
                self.setFormat(i, j - i, STYLES['keyword'])

                com = text[i:j]
                args_i = hit['args']
                i = args_i[0]

                # arguments
                for aa, j in zip(cmd.auto_arg, args_i[1:]):
                    aa = aa.get(com)
                    if aa is not None:
                        keywords = aa[0]().keywords
                        for m in re.finditer(r'\S+', text[i:j]):
                            mi = i + m.start()
                            mj = i + m.end()
                            arg = text[mi:mj]
                            if arg in keywords:
                                self.setFormat(mi, mj - mi, STYLES['argument'])

                    i = j + 1

                # warn about semantic errors
                if com in ('@', 'run'):
                    mi, mj = args_i[:2]
                    filename = text[mi:mj]
                    badext = '.py' if com == '@' else '.pml'
                    if filename.rstrip().endswith(badext):
                        self.setFormat(mi, mj - mi, STYLES['error'])

                self.setCurrentBlockState(multiline_state.get(com, -1))

        if text.endswith('\\'):
            if state_python != -1:
                self.setCurrentBlockState(CONTINUED_PYTHON | state_python)
            else:
                self.setCurrentBlockState(CONTINUED)
