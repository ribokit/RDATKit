# RDATKit (RNA Dataset ToolKIT)

**RDATkit** is a set of tools for parsing, analyzing, and publishing data of RNA chemical footprinting assays. It allows researchers to share their data using community standard formats, and helps them publish their results on indexable and shareable databases.

**RDATKit** is a package provides a set of *Python* and *MATLAB* scripts that facilitate saving and loading data to and from files with **RDAT** format. It also supports the **ISATAB** file format.


## Installation

#### MATLAB

- Download the zip or tar file of the repository and unpack; or 
```bash
git clone https://github.com/hitrace/RDATKit.git
```

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/RDATKit/MATLAB/`.

#### Python

To install **RDATKit**, simply:

- Copy `path.py.example` into `rdatkit/path.py`. Edit `rdatkit/path.py` following the instructions in the file to point to local installations of `RNAstructure`, `ViennaRNA`, and `VARNA`.

- Run:
```bash
cd path/to/RDATKit/
python setup.py install
```

For system-wide installation, you must have permissions and use with `sudo`.

**RDATKit** requires the following *Python* packages as dependencies, all of which can be installed through [`pip`](https://pip.pypa.io/).
```json
numpy >= 1.8.0
scipy >= 0.13.0
xlrd >= 0.9.2
xlwt >= 1.0.0
```

## Documentation

Documentation is available at https://hitrace.github.io/RDATKit/.

## License

Copyright &copy; of **RDATKit** _Source Code_ is described in [LICENSE.md](https://github.com/hitrace/RDATKit/blob/master/LICENSE.md).

<br/>
Developed by **Das lab**, _Leland Stanford Junior University_.
<br/>
README by [**t47**](http://t47.io/), *March 2016*.
