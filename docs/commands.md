# Commands added by pymolish-studio

Commands that PyMOL does not have. Each one works on the command line, in
scripts as `cmd.<name>`, and has `help <name>`. Since version 4.11.0 they
replace the commands of the pymolish plugin; remove the plugin from
`~/.pymol/startup/`, or its commands of the same name win.

## Groups

| Command | What it does | Example |
| --- | --- | --- |
| `group_next [group]` | Shows the next member of a group and hides the others. Key: Ctrl-N. | `group_next poses` |
| `group_previous [group]` | The same, backward. Key: Ctrl-P. | `group_previous` |
| `group_enable group [, each]` | Shows the members of a group, or every n-th one. | `group_enable poses, 10` |
| `group_disable group [, each]` | Hides the members; the group stays enabled. | `group_disable poses` |
| `copy_group group [, new_group [, suffix]]` | Copies a group with all members, groups inside it included. | `copy_group poses, backup, _bak` |
| `sort_group group [, reverse]` | Sorts the members by name, ignoring case. | `sort_group poses` |
| `extract_chain selection [, chains [, group]]` | Makes one new object per chain, `<object>_<chain>`. The source is not changed. | `extract_chain poses, A+B, poses_AB` |

The group row in the Objects panel has four of them in its A menu:
duplicate group, sort members (a second click sorts Z to A), enable members,
disable members.

## Many files at once

| Command | What it does | Example |
| --- | --- | --- |
| `load_files [directory [, format [, group [, pattern [, recursive [, prefix]]]]]]` | Loads every file of one format from a folder, one object per file. Unreadable files are reported and skipped; no object is replaced. | `load_files ~/poses, sdf, poses` |
| `export_objects selection [, directory [, format [, overwrite]]]` | Saves every object of a selection (a group, a name pattern, any selection) to `<object>.<format>`. Files that exist are skipped unless `overwrite=1`. | `export_objects poses, ~/out, pdb` |

In the window: File > Open Folder... and File > Export Objects.... Both show
the command they will run.

## Other

| Command | What it does | Example |
| --- | --- | --- |
| `palette [name [, recolor]]` | Sets the colors that new objects cycle through. | `palette vivid` |
