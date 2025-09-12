import sys
import time
from collections import defaultdict
from itertools import chain
from os.path import dirname
from pathlib import Path
from subprocess import PIPE, Popen


def create_all(generated_dir: str, pymol_dir: str = "."):
    """
    Generate various stuff
    """
    generated_dir_path = Path(generated_dir)
    pymol_dir_path = Path(pymol_dir)

    generated_dir_path.mkdir(parents=True, exist_ok=True)
    pymol_dir_path.mkdir(parents=True, exist_ok=True)

    create_shadertext(
        shader_dir=generated_dir_path / "data" / "shaders",
        shader_dir2=pymol_dir_path,
        output_header=generated_dir_path / "ShaderText.h",
        output_source=generated_dir_path / "ShaderText.cpp",
    )
    create_buildinfo(generated_dir, pymol_dir)


def create_shadertext(
    shader_dir: Path,
    shader_dir2: Path,
    output_header: Path,
    output_source: Path,
):
    varname = "_shader_cache_raw"
    include_deps = defaultdict(set)
    ifdef_deps = defaultdict(set)
    extension_regexp = {".gs", ".vs", ".fs", ".shared", ".tsc", ".tse"}

    # get all *.gs *.vs *.fs *.shared from the two input directories
    shader_files = set(
        path
        for path in chain(shader_dir.glob("**/*"), shader_dir2.glob("**/*"))
        if path.suffix in extension_regexp
    )

    with (
        open(output_header, "w") as output_header_file,
        open(output_source, "w") as output_source_file,
    ):
        output_header_file.write(f"extern const char * {varname}[];\n")
        output_source_file.write(f"const char * {varname}[] = {{\n")

        for shader_file in shader_files:
            output_source_file.write(f'"{shader_file.name}", ""\n')

            contents = shader_file.read_text()
            for line in contents.splitlines():
                line = line.strip()

                # skip blank lines and obvious comments
                if not line or line.startswith("//") and not "*/" in line:
                    continue

                # write line, quoted, escaped and with a line feed
                escaped_line = line.replace("\\", "\\\\").replace('"', r"\"")
                output_source_file.write(f'"{escaped_line}\\n"\n')

                # include and ifdef dependencies
                if line.startswith("#include"):
                    include_deps[line.split()[1]].add(shader_file.name)
                elif line.startswith("#ifdef") or line.startswith("#ifndef"):
                    ifdef_deps[line.split()[1]].add(shader_file.name)

            output_source_file.write(",\n")
        output_source_file.write("0};\n")

        # include and ifdef dependencies
        for varname, deps in [
            ("_include_deps", include_deps),
            ("_ifdef_deps", ifdef_deps),
        ]:
            output_header_file.write(f"extern const char * {varname}[];\n")
            output_source_file.write(f"const char * {varname}[] = {{\n")
            for name, item_deps in deps.items():
                item_deps = '", "'.join(sorted(item_deps))
                output_source_file.write(f'"{name}", "{item_deps}", 0,\n')
            output_source_file.write("0};\n")


def create_buildinfo(output_dir_path: str, pymoldir: str = "."):
    output_dir = Path(output_dir_path)
    sha_raw = Popen(["git", "rev-parse", "HEAD"], cwd=pymoldir, stdout=PIPE).stdout
    sha = sha_raw.read().strip().decode() if sha_raw is not None else ""

    info_file = output_dir / "PyMOLBuildInfo.h"
    info_file.write_text(
        f"""
    #define _PyMOL_BUILD_DATE {time.time()}
    #define _PYMOL_BUILD_GIT_SHA "{sha}"
    """
    )


if __name__ == "__main__":
    create_shadertext(
        shader_dir=Path(sys.argv[1]),
        shader_dir2=Path(sys.argv[2]),
        output_header=Path(sys.argv[3]),
        output_source=Path(sys.argv[4]),
    )
    create_buildinfo(dirname(sys.argv[4]), dirname(dirname(sys.argv[1])))
