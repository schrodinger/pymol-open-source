'''
Lighting Settings Plugin

(c) 2013-2017 Schrodinger Inc.
'''


def lightingsettings():
    from . import main
    return main.lightingsettings()


def __init_plugin__(self=None):
    from pymol import plugins
    plugins.addmenuitemqt('Lighting Settings', lightingsettings)
