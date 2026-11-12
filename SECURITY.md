# Security Policy

## Reporting a Vulnerability

Please report security vulnerabilities by email to
**[security+pymol@schrodinger.com](mailto:security+pymol@schrodinger.com)**.

Do not report security vulnerabilities through public GitHub issues, pull
requests, or discussions. Email is the only supported channel for security
reports.

Please include:

- The affected PyMOL version and platform (the `version` command in PyMOL
  prints the version).
- Steps to reproduce the issue. If the issue is triggered by a file (for
  example `.pse`, `.pdb`, or `.cif`), attach a sample file.
- The impact you expect, such as code execution, a crash, or information
  disclosure.

## What to Expect

The Schrödinger Information Security team receives and triages reports sent to
this address. We will acknowledge your report, work with the PyMOL developers
to confirm and fix the issue, and coordinate the disclosure timeline with you.
We credit reporters in the advisory unless you ask us not to.

Schrödinger does not run a bug bounty program for PyMOL and does not offer
monetary rewards for vulnerability reports.

## Supported Versions

Security fixes are made on the `master` branch and included in the next
release. Please check that the issue reproduces on the
[latest release](https://github.com/schrodinger/pymol-open-source/releases) or
on `master` before reporting.

## Scope

PyMOL is a scriptable application, and some features run code by design. The
following are not considered vulnerabilities:

- Running a Python (`.py`) or PyMOL (`.pml`) script that the user chose to run.
- Installing or loading a plugin.
- Using the `run`, `python`, or similar commands to execute code.

Anything that runs code or accesses data without the user asking for it, for
example when loading a structure, map, or session file, is in scope.
