# A* -------------------------------------------------------------------
# B* This file contains source code for the PyMOL computer program
# C* Copyright (c) Schrodinger, LLC.
# D* -------------------------------------------------------------------
# E* It is unlawful to modify or remove this copyright notice.
# F* -------------------------------------------------------------------
# G* Please see the accompanying LICENSE file for further information.
# H* -------------------------------------------------------------------
# I* Additional authors of this source file include:
# -*
# -*
# -*
# Z* -------------------------------------------------------------------

from typing import Iterable, Optional
from collections import defaultdict
from pymol import parsing


class Shortcut:
    def __init__(
        self,
        keywords: Optional[Iterable] = None,
        filter_leading_underscore: bool = True,
    ):
        keywords = list(keywords) if keywords is not None else []
        self.filter_leading_underscore = filter_leading_underscore
        self.keywords = (
            [keyword for keyword in keywords if keyword[:1] != "_"]
            if filter_leading_underscore
            else keywords
        )
        self.shortcut: dict[str, str | int] = {}
        self.abbreviation_dict = defaultdict(list)

        for keyword in self.keywords:
            self.optimize_symbols(keyword)

        self._rebuild_finalize()

    def __contains__(self, keyword: str) -> bool:
        return keyword in self.shortcut

    def __getitem__(self, keyword: str) -> Optional[int | str]:
        return self.shortcut.get(keyword)

    def __delitem__(self, keyword: str) -> None:
        self.keywords.remove(keyword)
        self.rebuild()

    def make_abbreviation(self, s: str, groups_length: int) -> str:
        """
        Example 1:
        Input: s:'abc_def_ghig', groups_length: 1
        Output: 'a_d_ghig'
        Example 2:
        Input: s:'abc_def', groups_length: 2
        Output: 'a_def'
        """
        groups = s.split("_")
        groups[:-1] = [c[0:groups_length] for c in groups[:-1]]
        return "_".join(groups)

    def optimize_symbols(self, keyword: str) -> None:
        for i in range(1, len(keyword)):
            substr = keyword[0:i]
            self.shortcut[substr] = 0 if substr in self.shortcut else keyword

        if "_" not in keyword:
            return

        for n in (1, 2):
            abbreviation = self.make_abbreviation(keyword, n)

            if keyword == abbreviation:
                continue

            self.abbreviation_dict[abbreviation].append(keyword)

            for i in range(abbreviation.find("_") + 1, len(abbreviation)):
                sub = abbreviation[0:i]
                self.shortcut[sub] = 0 if sub in self.shortcut else keyword

    def rebuild(self, keywords: Optional[Iterable] = None) -> None:
        keywords = list(keywords) if keywords is not None else []
        self.keywords = (
            [keyword for keyword in keywords if keyword[:1] != "_"]
            if self.filter_leading_underscore
            else keywords
        )
        # optimize symbols
        self.shortcut = {}
        self.abbreviation_dict = defaultdict(list)
        for keyword in self.keywords:
            self.optimize_symbols(keyword)

        self._rebuild_finalize()

    def _rebuild_finalize(self) -> None:
        for abbreviation, keywords in self.abbreviation_dict.items():
            if len(keywords) == 1:
                self.shortcut[abbreviation] = keywords[0]
        for keyword in self.keywords:
            self.shortcut[keyword] = keyword

    def interpret(
        self, keyword: str, mode: bool = False
    ) -> Optional[int | str | list[str]]:
        """
        Returns None (no hit), str (one hit) or list (multiple hits)

        keyword = str: query string, setting prefix or shortcut
        mode = True/False: if mode=True, do prefix search even if kee has exact match
        """
        if keyword == "":
            return self.keywords

        result = self.shortcut.get(keyword)
        if result is None:
            return
        if result and not mode:
            return result

        # prefix search
        unique_keywords = set(
            word for word in self.keywords if word.startswith(keyword)
        )
        for abbreviation, keywords in self.abbreviation_dict.items():
            if abbreviation.startswith(keyword):
                unique_keywords.update(keywords)
        # no match
        if not unique_keywords:
            return

        # single match: str
        # multiple matches: list
        return (
            unique_keywords.pop()
            if len(unique_keywords) == 1
            else list(unique_keywords)
        )

    def append(self, keyword) -> None:
        self.keywords.append(keyword)
        self.optimize_symbols(keyword)
        self._rebuild_finalize()

    def auto_err(
        self, keyword: str, descrip: Optional[str] = None
    ) -> Optional[int | str | list[str]]:
        if keyword == "":
            return

        result = self.interpret(keyword)

        if result is None and descrip is not None:
            msg = f"Error: unknown {descrip}: '{keyword}'."
            lst = self.interpret("")
            if isinstance(lst, list) and len(lst) < 100:
                lst.sort()
                lst = parsing.list_to_str_list(lst)
                msg += " Choices:\n" + "\n".join(lst)
                raise parsing.QuietException(msg)

        if isinstance(result, list) and descrip is not None:
            lst = parsing.list_to_str_list(result)
            options = "\n".join(lst)
            msg = f"Error: ambiguous {descrip}\\n {options}"
            raise parsing.QuietException(msg)

        return result
