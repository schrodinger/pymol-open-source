# Adopted from pillow and pycapnp

import sys

from setuptools.build_meta import build_wheel

backend_class = build_wheel.__self__.__class__


class _CustomBuildMetaBackend(backend_class):
    def run_setup(self, setup_script="setup.py"):
        if self.config_settings:
            flags = [
                f"--{k}={v}"
                for k, v in self.config_settings.items()
            ]
            sys.argv = sys.argv[:1] + ["build_ext"] + flags + sys.argv[1:]
        return super().run_setup(setup_script)

    def build_wheel(
        self, wheel_directory, config_settings=None, metadata_directory=None
    ):
        self.config_settings = config_settings
        return super().build_wheel(wheel_directory, config_settings, metadata_directory)

    def build_editable(
        self, wheel_directory, config_settings=None, metadata_directory=None
    ):
        self.config_settings = config_settings
        return super().build_editable(
            wheel_directory, config_settings, metadata_directory
        )

_backend = _CustomBuildMetaBackend()
build_wheel = _backend.build_wheel
build_editable = _backend.build_editable
