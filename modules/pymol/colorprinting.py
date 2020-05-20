##
## (c) Schrodinger, Inc.
##

import re


def htmlspecialchars(s):
    '''Convert special characters to HTML entities
    '''
    return (s
        .replace('&', '&amp;')
        .replace('<', '&lt;')
        .replace('>', '&gt;')
    )

def text2html(s, whitespace=True):
    '''Convert plain text to HTML
    '''
    s = htmlspecialchars(s)

    if whitespace:
        s = s.replace(' ', '&nbsp;').replace('\n', '<br>')

    return s


error = print
warning = print
suggest = print
parrot = print


def print_exc(strip_filenames=()):
    '''Print a colored traceback, with an (optionally) reduced stack
    '''
    import sys, traceback

    exc_type, exc_value, tb = sys.exc_info()

    for filename in strip_filenames:
        while tb and filename.startswith(tb.tb_frame.f_code.co_filename):
            tb = tb.tb_next

    s_list = traceback.format_exception(exc_type, exc_value, tb)
    error(''.join(s_list).rstrip())

    return exc_type, exc_value, tb
