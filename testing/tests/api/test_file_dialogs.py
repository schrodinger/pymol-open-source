"""
Tests for the Qt file dialogs (modules/pmg_qt/file_dialogs.py).
"""
from unittest import mock

import pytest

pytest.importorskip('pymol.Qt')

from pmg_qt import file_dialogs


def _run_movie_export(fmt, ray, fname):
    """
    Open the "Export Movie" dialog with a mocked parent and form, select
    `fmt` and the ray/draw option, and click "OK". Returns the mocked parent.
    """
    parent = mock.Mock()
    parent.cmd.get_viewport.return_value = (640, 480)
    parent.cmd.get_setting_int.return_value = 90

    form = parent.load_form.return_value
    form.input_encoder.count.return_value = 2
    form.input_encoder.itemText.return_value = 'ffmpeg'
    form.input_encoder.currentText.return_value = 'ffmpeg'
    form.input_ray.isChecked.return_value = ray
    form.input_width.value.return_value = 640
    form.input_height.value.return_value = 480
    form.input_quality.value.return_value = 90
    for name in ('png', 'mp4', 'mpg', 'mov', 'gif'):
        getattr(form, 'format_' + name).isChecked.return_value = (name == fmt)

    with mock.patch('pymol.movie.find_exe', return_value='ffmpeg'):
        file_dialogs.file_save_mpeg(parent, 'png' if fmt == 'png' else None)

    run = form.button_ok.clicked.connect.call_args[0][0]

    with mock.patch.object(file_dialogs, 'getSaveFileNameWithExt',
                           return_value=fname):
        run()

    return parent


@pytest.mark.parametrize('ray, mode', [(True, 2), (False, 1)])
def test_movie_export_png_is_modal(ray, mode, tmp_path):
    # GitHub #519: modal=-1 made ray-traced export block the GUI event loop
    parent = _run_movie_export('png', ray, str(tmp_path / 'frame.png'))

    parent.cmd.mpng.assert_called_once()
    kwargs = parent.cmd.mpng.call_args.kwargs
    assert kwargs['mode'] == mode
    assert kwargs['modal'] == 1


@pytest.mark.parametrize('ray, mode', [(True, 'ray'), (False, 'draw')])
def test_movie_export_video_is_modal(ray, mode, tmp_path):
    # GitHub #519: same freeze for MP4/GIF/MOV export through movie.produce
    parent = _run_movie_export('mp4', ray, str(tmp_path / 'movie.mp4'))

    parent.cmd.movie.produce.assert_called_once()
    kwargs = parent.cmd.movie.produce.call_args.kwargs
    assert kwargs['mode'] == mode
    assert kwargs['modal'] == 1
