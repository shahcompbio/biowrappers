[toc]

# Overview

The `biowrappers` package is a pure Python package, so it can be installed easily.
However, `biowrappers` requires a diverse set of tools to run, some of which have special wrapper scripts to work with the library.
To install these tools we will leverage the `conda` package manager.

# Install from binary

TODO: Find a way to upload conda packages somewhere accessible.

# Build conda binaries

## Install conda

Download [`conda`](http://conda.pydata.org/miniconda.html) and follow install directions.
Note the prefix where you install for later.

## Install `conda-build`

```
conda install conda-build
```

## Build `conda` packages

```
git clone git@bitbucket.org:aroth85/biowrappers.git biowrappers

cd biowrappers/recipes

conda build biowrappers-stack
```

Built packages will be in `CONDA_INSTALL_DIR/conda-bld/PLATFORM`. 
Where CONDA_INSTALL_DIR is the path to PREFIX you used when installing `conda` and PLATFORM is you platform.

## Install packages

```
conda install biowrappers-stack -c file://my/conda/path/conda-bld
```

Replace /my/conda/path with the absolute path of CONDA_INSTALL_DIR.