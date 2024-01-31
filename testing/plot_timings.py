"""
Create plots for timings

Will create a temporary directory and put plots in PNG format there,
together with an index.html file.

"""

from __future__ import print_function

import os, sys, tempfile, re, socket, time
import datetime
from optparse import OptionParser
from collections import defaultdict
from matplotlib import pyplot, rcParams, dates

# plot setup
rcParams['figure.figsize'] = 5.0, 2.5
rcParams['font.size'] = 9

# command line options
parser = OptionParser()
parser.add_option("-b", "--browse", action="store_true", dest="browse")
parser.add_option("-q", "--quiet", action="store_true", dest="quiet")
options = parser.parse_args()[0]

# file with timing results
tabname = os.getenv("PYMOLTESTTIMINGS", "timings.tab")

# read file
db = defaultdict(lambda: defaultdict(list))
for line in open(tabname):
    a = line.rstrip('\n').split("\t")
    timestamp = float(a[0])
    try:
        mac = a[9] + a[5]
    except:
        mac = a[1]
        mac = ':'.join(mac[i:i+2] for i in range(0, len(mac), 2))
    value = float(a[2])
    key = '%s(%s)' % (a[3], a[4])
    db[key][mac].append((timestamp, value))

# helper function for unique PNG filenames
used_png = set()
def get_unused_png(key):
    r = key = re.sub(r'[^-\w.]', '_', key)
    i = 0
    while r in used_png:
        i += 1
        r = key + '-%d' % i
    return r + '.png'

# create output dir
outdir = tempfile.mkdtemp()
htmlout = open(os.path.join(outdir, "index.html"), "w")
print("<h1>PyMOL Benchmarks,", socket.gethostname(), file=htmlout)
print(time.strftime("%D-%T"), "</h1>", file=htmlout)

# make plots
for key in sorted(db):
    data = db[key]
    pyplot.clf()
    fig, ax = pyplot.subplots(1)
    maxy = 0.0
    for mac in data:
        x, y = zip(*data[mac])
        x = list(map(datetime.datetime.fromtimestamp, x))
        maxy = max(maxy, max(y))
        ax.plot(x, y, "o-", label=mac[:290])

    fig.autofmt_xdate()
    ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%m-%d'))
    ax.yaxis.set_label_text("Seconds")
    ax.set_title(key)

    pyplot.grid(True)
    pyplot.ylim(0, maxy * 1.1)
    pyplot.legend(loc="best")
    pngname = get_unused_png(key)
    pyplot.savefig(os.path.join(outdir, pngname), dpi=70)
    print("<img src='%s'>" % (pngname), file=htmlout)

# done
if not options.quiet:
    print(outdir)

# open index.html
if options.browse:
    outhtml = os.path.join(outdir, "index.html")
    if sys.platform.startswith("darwin"):
        os.system("open " + outhtml)
    else:
        import webbrowser
        webbrowser.open(outhtml)
