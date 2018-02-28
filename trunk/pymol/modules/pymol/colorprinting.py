##
## (c) Schrodinger, Inc.
##

from __future__ import absolute_import as _
from __future__ import print_function as _
from __future__ import division as _

import re

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
