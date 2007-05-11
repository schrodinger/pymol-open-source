from compileall import compile_dir

compile_dir("modules",force=True)
compile_dir("epymol/modules",force=True)
compile_dir("ipymol/modules",force=True)

