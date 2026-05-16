# RDATKit (RNA Dataset ToolKIT)

**RDATkit** is a set of tools for parsing, analyzing, and publishing data of RNA chemical footprinting assays. It allows researchers to share their data using community standard formats, and helps them publish their results on indexable and shareable databases.

**RDATKit** is a package provides a set of *Python* and *MATLAB* scripts that facilitate saving and loading data to and from files with **RDAT** format. It also supports the **ISATAB** file format.


## Installation

#### MATLAB

- Download the zip or tar file of the repository and unpack; or 
```bash
git clone https://github.com/ribokit/RDATKit.git
```

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/RDATKit/MATLAB/`.

#### Python

To install **RDATKit**, simply:

```bash
pip install rdat_kit
```

**RDATKit** requires the following *Python* packages as dependencies, all of which can be installed through [`pip`](https://pip.pypa.io/).
```json
numpy >= 1.8.0
scipy >= 0.13.0
xlrd >= 0.9.2
xlwt >= 1.0.0
```

## Command-line interface

As of v1.7.0, `rdat_kit` registers a CLI with two subcommands:

```bash
# Validate one or more RDAT files (parses + runs RDATFile.validate())
rdat_kit validate path/to/entry.rdat

# Emit a Jekyll/RMDB front-matter .md stub for use in
# https://github.com/DasLab/rmdb.github.io
rdat_kit to_md path/to/entry.rdat > _entries/RMDB_ID.md
```

`validate` exits with code 0 on success, 1 on parse failure, 2 if
RDATFile.validate() returned warnings. `to_md` derives the RMDB_ID
from the filename (`<PREFIX>_<CHEM>_<NNNN>` pattern) and accepts
`--rmdb-id` to override.

## Documentation

Documentation is available at https://ribokit.github.io/RDATKit/.

## License

Copyright &copy; of **RDATKit** _Source Code_ is described in [LICENSE.md](https://github.com/ribokit/RDATKit/blob/master/LICENSE.md).

<hr/>

Developed by [**Das lab**](https://daslab.stanford.edu), _Leland Stanford Junior University_.