import pytest

from pymol.shortcut import Shortcut


@pytest.fixture
def sc() -> Shortcut:
    return Shortcut(["foo", "bar", "baz", "com", "com_bla", "com_xxx"])


@pytest.mark.parametrize(
    "keyword, expected_result",
    [
        ("a", False),
        ("w", True),
        ("war", True),
    ],
)
def test_contains(keyword: str, expected_result: bool):
    shortcut = Shortcut(["warren", "wasteland", "electric", "well"])
    assert (keyword in shortcut) is expected_result


def test_interpret():
    shortcut = Shortcut(["warren", "wasteland", "electric", "well"])
    list_result = shortcut.interpret("w")
    assert list_result is not None
    assert not isinstance(list_result, int)
    assert sorted(list_result) == ["warren", "wasteland", "well"]

    string_result = shortcut.interpret("e")
    assert list_result is not None
    assert string_result == "electric"


def test_all_keywords(sc: Shortcut):
    assert ["foo", "bar", "baz", "com", "com_bla", "com_xxx"] == sc.interpret("")


@pytest.mark.parametrize(
    "prefixs, expected_result",
    [
        (["f", "fo", "foo"], "foo"),
        (["b", "ba"], ["bar", "baz"]),
        (["bar"], "bar"),
        (["c", "co"], ["com", "com_bla", "com_xxx"]),
        (["com"], "com"),
    ],
)
def test_full_prefix_hits(
    sc: Shortcut, prefixs: list[str], expected_result: str | list[str]
):
    for prefix in prefixs:
        result = sc.interpret(prefix)
        result = sorted(result) if isinstance(result, list) else result
        assert expected_result == result


def test_append(sc: Shortcut):
    sc.append("foo_new")

    assert ["foo", "foo_new"], sc.interpret("f")
    assert "foo", sc.interpret("foo")
    assert "foo_new", sc.interpret("foo_")

    assert "" not in sc


def test_abbreviations(sc: Shortcut):
    sc.append("foo_new")

    assert "foo_new" == sc.interpret("f_")
    assert "foo_new" == sc.interpret("f_new")
    assert "foo_new" == sc.interpret("fo_")
    assert "com_xxx" == sc.interpret("c_x")
    assert "com_xxx" == sc.interpret("c_xxx")
    assert "com_xxx" == sc.interpret("co_x")


def test_missing_key(sc: Shortcut):
    assert None is sc.interpret("missing_key")


def test_auto_error(sc: Shortcut):
    assert None is sc.auto_err("")
    assert None is sc.auto_err("missing_key")

    result = sc.auto_err("co")
    assert isinstance(result, list)
    assert ["com", "com_bla", "com_xxx"] == sorted(result)
    assert "com", sc.auto_err("com")


def test_interpret_mode_true(sc: Shortcut):
    assert "foo" == sc.interpret("f", True)

    result = sc.interpret("com", True)
    assert isinstance(result, list)
    assert ["com", "com_bla", "com_xxx"] == sorted(result)

    sc.append("foo_new")
    result = sc.interpret("foo", True)
    assert isinstance(result, list)
    assert ["foo", "foo_new"] == sorted(result)


def test_rebuild(sc: Shortcut):
    coms = ["com", "com_bla", "com_xxx"]
    sc.rebuild(coms)

    assert None is sc.interpret("f")
    assert None is sc.interpret("foo")

    result = sc.interpret("c")
    assert isinstance(result, list)
    assert coms == sorted(result)

    result = sc.interpret("com", True)
    assert isinstance(result, list)
    assert coms == sorted(result)

    assert "com" == sc.interpret("com")
    assert "com_xxx" == sc.interpret("c_x")
