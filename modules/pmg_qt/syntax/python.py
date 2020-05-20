import re
import keyword

from . import *

STYLES = {
    'keyword': textformat('darkOrange'),
    'builtins': textformat('green'),
    'decorator': textformat('darkBlue'),
    'defclass': textformat('darkGreen'),
    'string': textformat('magenta'),
    'comment': textformat('blue'),
    'self': textformat('black', 'bold'),
    'numbers': textformat('brown'),
}


# enum for block state
class QUOTESTATE:
    NONE = 0
    DOUBLE = 1
    SINGLE = 2
    TRIPLEOFFSET = 2


# quotes (QUOTESTATE indices)
QUOTES = [None, '"', "'", '"""', "'''"]


def find_first(subs, text, start):
    '''Index of first occurence of `subs` in `text`. Return None if not found.'''
    best = None
    for sub in subs:
        index = text.find(sub, start, best)
        if index != -1:
            best = index
    return best


python_rules = [
    (re.compile(pat), index, fmt)
    for (pat, index, fmt) in [
        # decorators
        (r'^\s*@[\w\.]+', 0, STYLES['decorator']),

        # built in names
        (r'(?<!\.)\b(?:' + '|'.join(__builtins__) + r')\b', 0, STYLES['builtins']),

        # 'self'
        (r'\bself\b', 0, STYLES['self']),

        # 'def' or 'class' followed by an identifier
        (r'\b(?:def|class)\s+(\w+)', 1, STYLES['defclass']),

        # Numeric literals
        (r'\b[+-]?(?:[0-9]+|0[xX][0-9A-Fa-f]+|0[bB][01]+)[lL]?\b', 0, STYLES['numbers']),
        (r'\b[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\b', 0, STYLES['numbers']),

        # keywords
        (r'\b(?:' + '|'.join(keyword.kwlist) + r')\b', 0, STYLES['keyword']),
    ]
]

re_end_of_line = re.compile(r'$', re.M)


def highlightPython(self, state, text, offset=0, end=None):
    '''
    @type self: QSyntaxHighlighter
    @type state: int
    @type text: str
    @type offset: int
    @type end: int or None

    @rtype: int
    @return current block state
    '''
    # support highlighting substrings
    text = text[offset:end]

    # regular expressions
    for expression, nth, format in python_rules:
        for match in expression.finditer(text):
            start, index = match.span(nth)
            self.setFormat(start + offset, index - start, format)

    start = 0
    index = 0
    L = len(text)

    while index < L:
        if state == QUOTESTATE.NONE:
            # find opening quote
            start = find_first(('"', "'", '#'), text, index)

            if start is None:
                # nothing found
                break

            index = start + 1

            if text[start] == '#':
                # found comment
                index += re_end_of_line.search(text[index:]).end()
                self.setFormat(start + offset, index - start,
                               STYLES['comment'])
                continue

            if text[start] == '"':
                # found double quote
                state = QUOTESTATE.DOUBLE
            else:
                # found single quote
                state = QUOTESTATE.SINGLE

            # check for triple quote
            if text[index:index + 2] == QUOTES[state] * 2:
                state += QUOTESTATE.TRIPLEOFFSET
                index += 2

        quote = QUOTES[state]
        quotelen = len(quote)

        # find closing quote
        while index < L:
            if text[index] == '\\':
                index += 1
            elif text[index:index + quotelen] == quote:
                index += quotelen
                state = QUOTESTATE.NONE
                break
            index += 1

        # raw/unicode/bytes prefix
        if start > 0 and text[start - 1].lower() in 'rub':
            start -= 1

        self.setFormat(start + offset, index - start, STYLES['string'])

    return state


class Highlighter(QtGui.QSyntaxHighlighter):
    def highlightBlock(self, text):
        state = max(0, self.previousBlockState())
        state = highlightPython(self, state, text)
        self.setCurrentBlockState(state)
