"""
README and license files (plan step 10): relative links resolve, LICENSE
still ends with the upstream notice, and the PyMOL(TM) slogans are used the
way LICENSE allows.
"""

import hashlib
import os
import re

from pymol import testing

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
DOCS = ["README.md", "CONTRIBUTING.md"] + sorted(
    "docs/" + name
    for name in os.listdir(os.path.join(ROOT, "docs"))
    if name.endswith(".md")
)
# every place that shows a slogan to users
SLOGAN_FILES = DOCS + ["pyproject.toml", ".github/workflows/build.yml"]
SLOGANS = (
    "Includes PyMOL(TM).",
    "Based on PyMOL(TM) technology.",
    "Contains PyMOL(TM) source code.",
    "Built using PyMOL(TM).",
)
FOOTNOTE = "PyMOL is a trademark of Schrodinger, LLC."
UPSTREAM_LENGTH = 2268
UPSTREAM_SHA256 = "a1fd32ff4f11d3b7eedd8dd4ecfcec1bd65374e16ee8a823899dc3c2be13d85e"


def read(name):
    # text mode: a Windows checkout may have CRLF line ends
    with open(os.path.join(ROOT, name), encoding="utf-8") as handle:
        return handle.read()


class TestLicenseDocs(testing.PyMOLTestCase):
    def test_links(self):
        for name in DOCS:
            targets = re.findall(r"\]\(([^)#\s]+)", read(name))
            targets += re.findall(r'src="([^"]+)"', read(name))
            for target in targets:
                if "://" in target:
                    continue
                path = os.path.join(ROOT, os.path.dirname(name), target)
                self.assertTrue(os.path.exists(path), "%s: %s" % (name, target))

    def test_upstream_notice(self):
        tail = read("LICENSE")[-UPSTREAM_LENGTH:]
        self.assertEqual(hashlib.sha256(tail.encode()).hexdigest(), UPSTREAM_SHA256)
        self.assertIn("MIT License", read("LICENSE-MIT"))

    def test_slogans(self):
        for name in SLOGAN_FILES:
            text = read(name)
            self.assertNotRegex(text, r"(?i)pymolish(-studio)?\(TM\)", name)
            if "PyMOL(TM)" in text:
                self.assertIn(FOOTNOTE, text, name)
            for slogan in SLOGANS:
                if slogan in text:
                    self.assertEqual(
                        text.count(slogan),
                        text.count("pymolish-studio. " + slogan),
                        name,
                    )
        self.assertTrue(read("README.md").rstrip().endswith(FOOTNOTE))
