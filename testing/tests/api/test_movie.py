from unittest import mock

import pytest

from pymol import movie


@pytest.mark.parametrize('kwargs, expected_modal', [
    ({}, -1),
    ({'modal': 1}, 1),
    ({'modal': 0}, 0),
])
def test_produce_forwards_modal(kwargs, expected_modal, tmp_path, monkeypatch):
    # GitHub #519: the GUI needs movie.produce to render modally
    monkeypatch.setattr(movie, 'find_exe', lambda exe: exe)
    monkeypatch.setattr(movie, '_encode', mock.Mock())

    _self = mock.Mock()
    _self.get_modal_draw.return_value = 0
    _self.count_frames.return_value = 3

    movie.produce(str(tmp_path / 'movie.mp4'), mode='ray', width=64,
                  height=48, encoder='ffmpeg', quality=90, _self=_self,
                  **kwargs)

    _self.mpng.assert_called_once()
    assert _self.mpng.call_args.kwargs['modal'] == expected_modal
    assert _self.mpng.call_args.kwargs['mode'] == 2
