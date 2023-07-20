import pymol

pymol.__path__.append(".")
pymol.__path__.append("tests/helpers")

collect_ignore = [
    "tests/helpers",
]
