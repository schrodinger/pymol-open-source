# Contributing

pymolish-studio has one goal: **smaller, faster, simpler, prettier**. A change
is welcome when it serves one of those four words. Removing code counts as
much as adding it.

## Before you start

- **Small fix or removal:** open a pull request.
- **New feature, new dependency or a change to how a `cmd` command behaves:**
  open an issue first, so nobody builds something that will not be merged.
- The [roadmap](PLAN.md) shows what is planned; each step has a detailed plan
  in `docs/plans/`.

## Make the change

1. Build an editable install: [Install from source](docs/install.md#working-on-the-code).
2. Keep the change small: one topic per pull request.
3. Leave a test behind for anything that is not trivial, in `testing/tests/`
   (Python) or `layerCTest/` (C++).
4. Run the suite: `pymolish -ckqy testing/testing.py --run all`.
5. Format only the files you touched: `clang-format -i` for C++,
   `ruff format` for Python.
6. If the behavior of a `cmd` command changes, say so in the pull request
   and add a line to `ChangeLog`.

Commit subjects look like `fix(labels): atom labels were solid blocks`. The
type is one of `remove`, `perf`, `simplify`, `ui`, `fix`, `deps`, `ci`,
`docs`.

The full rules, written for people and for AI coding agents, are in
[AGENTS.md](AGENTS.md).

## License of your contribution

What you contribute is published under the [MIT license](LICENSE-MIT). You
keep your copyright; you do not transfer it to anyone.

Do not change or remove the copyright notices that are already in the files,
and do not edit the notice in [LICENSE](LICENSE). Why two licenses apply is
explained in [Licenses](docs/licenses.md).

The developer notice of the original project is kept, unchanged, in
[docs/DEVELOPERS-upstream.txt](docs/DEVELOPERS-upstream.txt). It describes
contributing to Schrodinger's project and does not apply here.
