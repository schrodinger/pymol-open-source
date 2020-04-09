#pragma once

/**
 * Copyable arguments container for ExecutiveLoad
 */
struct ExecutiveLoadArgs {
  std::string fname;
  std::string content;
  cLoadType_t content_format;
  std::string object_name;
  int state;
  int zoom;
  int discrete;
  int finish;
  int multiplex;
  int quiet;
  std::string plugin;
  std::string object_props;
  std::string atom_props;
  bool mimic;
  int plugin_mask = 0;
};

/**
 * Prepare arguments for ExecutiveLoad
 */
pymol::Result<ExecutiveLoadArgs> ExecutiveLoadPrepareArgs(PyMOLGlobals* G,
    pymol::null_safe_zstring_view fname, const char* content,
    int content_length, cLoadType_t content_format,
    const char* object_name_proposed, int state, int zoom, int discrete,
    int finish, int multiplex, int quiet, const char* plugin_arg,
    const char* object_props, const char* atom_props, bool mimic = false);

/**
 * Implementation of ExecutiveLoad
 *
 * @param args Prepared arguments
 * @return true on success
 */
pymol::Result<> ExecutiveLoad(PyMOLGlobals* G, ExecutiveLoadArgs const& args);
