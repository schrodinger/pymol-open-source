def process(sele):
    # convert object/index tuples into selection strings
    if isinstance(sele, tuple):
        return "%s`%d" % sele

    # convert unicode hyphens to dashes
    return str(sele).replace(u'\u2212', '-')
