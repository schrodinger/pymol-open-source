'''
Headless performance baseline. Run from the repository root:

    pymol -ckqy testing/bench.py

Prints a markdown table (median and min..max of N repeats, peak RSS so far)
for load, representation build + CPU ray, isosurface, and PSE round trip.
GUI frame latency is not measured here; it needs a real GL context.
'''
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time

from pymol import cmd

N = 5
# pymol execs scripts in its own namespace, so __file__ is not this file
DATA = os.path.join(os.path.dirname(os.path.abspath(sys.argv[-1])), 'data')
TMP = tempfile.mkdtemp(prefix='pymol-bench-')
rows = []


def bench(name, setup, run, check=None):
    times = []
    for _ in range(N):
        cmd.reinitialize()
        setup()
        t0 = time.perf_counter()
        run()
        times.append(time.perf_counter() - t0)
        if check:
            check()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    rows.append((name, statistics.median(times), min(times), max(times), rss))
    print(f'  {name}: {statistics.median(times):.3f} s', file=sys.stderr)


def load_1aon():
    cmd.load(os.path.join(DATA, '1aon.pdb.gz'), 'm')


def check_1aon():
    assert cmd.count_atoms('m') == 58870, cmd.count_atoms('m')


def load_map():
    cmd.load(os.path.join(DATA, 'emd_1155.ccp4'), 'map')


cmd.set('max_threads', os.cpu_count())

bench('load 1aon.pdb.gz (58870 atoms)', lambda: None, load_1aon, check_1aon)

bench('ray 800x600, cartoon (1aon)',
      lambda: (load_1aon(), cmd.hide('everything'), cmd.show('cartoon')),
      lambda: cmd.ray(800, 600))

bench('ray 800x600, surface (1aon)',
      lambda: (load_1aon(), cmd.hide('everything'), cmd.show('surface')),
      lambda: cmd.ray(800, 600))

# headless PyMOL builds the mesh lazily; a tiny ray forces it
bench('isosurface emd_1155.ccp4 level 2.0 (+ray 64x64)',
      load_map, lambda: (cmd.isosurface('iso', 'map', 2.0), cmd.ray(64, 64)))

bench('isosurface level change 2.0 -> 3.0 (+ray 64x64)',
      lambda: (load_map(), cmd.isosurface('iso', 'map', 2.0), cmd.ray(64, 64)),
      lambda: (cmd.isosurface('iso', 'map', 3.0), cmd.ray(64, 64)))

pse = os.path.join(TMP, 'bench.pse')
bench('save PSE (1aon, cartoon)',
      lambda: (load_1aon(), cmd.show('cartoon')),
      lambda: cmd.save(pse))
bench('load PSE (1aon, cartoon)',
      lambda: None, lambda: cmd.load(pse), check_1aon)

exe = os.path.join(os.path.dirname(sys.executable), 'pymol')
startup = []
for _ in range(N):
    t0 = time.perf_counter()
    subprocess.run([exe, '-ckq', '-d', 'quit'], check=True, stdout=subprocess.DEVNULL)
    startup.append(time.perf_counter() - t0)
rows.append(('startup `pymol -ckq -d quit`', statistics.median(startup),
             min(startup), max(startup), float('nan')))

print(f'\nCommit {subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()}, '
      f'PyMOL {cmd.get_version()[0]}, Python {platform.python_version()}, '
      f'{platform.system()} {platform.machine()}, {os.cpu_count()} CPUs, '
      f'PSE {os.path.getsize(pse) / 1e6:.1f} MB\n')
print('| Workload | median | min..max | peak RSS |')
print('| --- | --- | --- | --- |')
for name, med, lo, hi, rss in rows:
    print(f'| {name} | {med:.3f} s | {lo:.3f}..{hi:.3f} s | {rss:.0f} MB |')
