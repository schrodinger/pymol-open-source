<img src="data/pymol/icons/pymolish.png" height="100" align="right" alt="pymolish-studio icon" />

# pymolish-studio

pymolish-studio. Based on PyMOL(TM) technology.

A desktop program for looking at molecules in 3D: open a protein structure,
turn it around, color it, save a picture. It is
[Open-Source PyMOL](https://github.com/schrodinger/pymol-open-source) with a
cleaner window, less old code and faster heavy operations.

> pymolish-studio is an independent fork. It is not a Schrodinger product and
> is not affiliated with or endorsed by Schrodinger, LLC. The official PyMOL
> is at [pymol.org](https://pymol.org).

![The pymolish-studio window with a protein shown as cartoon](docs/img/window.png)

## Why use it

- **A window that is easy to read.** Dark theme, one toolbar with icons, a
  panel that lists your objects, and a size button for high-resolution
  screens.
- **Faster where it was slow.** A large session opens 6.6 times faster and a
  large surface builds in about half the time.
- **Nothing to relearn.** Commands, scripts, Qt plugins and files (PDB,
  mmCIF, SDF, PSE sessions) work as they do in PyMOL.
- **Less to carry.** Old interfaces and unused features are gone, so there is
  less to install and less to break.

The full list is in [What is different](docs/features.md).

## Try it

1. Download the archive for your system from the
   [Releases page](https://github.com/adrian-jasinski/pymolish-studio/releases).
   There is one for Linux, Windows and macOS (Apple silicon). You do not need
   Python.
2. Unpack it and start `pymolish` (`pymolish.exe` on Windows).
3. Choose File > Get PDB, type `1ubq`, and drag with the mouse.

[Getting started](docs/getting-started.md) takes it from there.

These are preview builds and they are not signed, so Windows and macOS may
warn before the first start. `SHA256SUMS` on the Releases page lists the
checksums.

## Documentation

| I want to | Read |
| --- | --- |
| Open a structure and make a picture | [Getting started](docs/getting-started.md) |
| Know what changed compared with PyMOL | [What is different](docs/features.md) |
| Build it from source | [Install from source](docs/install.md) |
| Write a plugin | [Writing a plugin](docs/plugins.md) |
| Help with the code | [Contributing](CONTRIBUTING.md) |
| See what is planned | [Roadmap](PLAN.md), [finished steps](PLAN-DONE.md), [ChangeLog](ChangeLog) |
| Know what I may do with it | [Licenses](docs/licenses.md) |

## License

pymolish-studio. Contains PyMOL(TM) source code.

Free to use, change and share, also at work. Two licenses apply, and a copy
of the software needs both files:

- Code from Open-Source PyMOL, Copyright (C) Schrodinger, LLC:
  [LICENSE](LICENSE).
- Code added by pymolish-studio: MIT, [LICENSE-MIT](LICENSE-MIT).

[Licenses](docs/licenses.md) explains this in plain words and lists the
third-party parts.

---

PyMOL is a trademark of Schrodinger, LLC.
