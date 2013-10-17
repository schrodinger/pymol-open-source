"""
Create plots for timings

Will create a temporary directory and put plots in PNG format there,
together with an index.html file.

"""

import os, sys, tempfile, re, socket, time
from optparse import OptionParser
from collections import defaultdict
from matplotlib import pyplot, rcParams

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
db = defaultdict(list)
for line in open(tabname):
    a = line.split("\t")
    timestamp = float(a[0])
    key = a[1]
    value = float(a[2])
    db[key].append((timestamp, value))

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
print >> htmlout, "<h1>PyMOL Benchmarks,", socket.gethostname()
print >> htmlout, time.strftime("%D-%T"), "</h1>"

# make plots
for key in sorted(db):
    x, y = zip(*db[key])
    pyplot.clf()
    pyplot.ylim(0, max(y) * 1.1)
    pyplot.title(key)
    pyplot.ylabel("Timing in Seconds")
    pyplot.grid(True)
    pyplot.plot(x, y)
    pngname = get_unused_png(key)
    pyplot.savefig(os.path.join(outdir, pngname), dpi=70)
    print >> htmlout, "<img src='%s'>" % (pngname)

# done
if not options.quiet:
    print outdir

# open index.html
if options.browse:
    outhtml = os.path.join(outdir, "index.html")
    if sys.platform.startswith("darwin"):
        os.system("open " + outhtml)
    else:
        import webbrowser
        webbrowser.open(outhtml)
