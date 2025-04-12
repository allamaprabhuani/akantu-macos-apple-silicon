# `Akantu`: Swiss-Made Open-Source Finite-Element Library

[![joss](https://joss.theoj.org/papers/3abf3c7945cc9a016a946ce9e02e357f/status.svg)](https://joss.theoj.org/papers/3abf3c7945cc9a016a946ce9e02e357f) [![license](https://img.shields.io/badge/license-LGPLv3-green)](https://www.gnu.org/licenses/lgpl-3.0.en.html) [![readthedoc](https://readthedocs.org/projects/akantu/badge/?version=master)](https://akantu.readthedocs.io/en/latest/?badge=master)

`Akantu` means a little element in Kinyarwanda, a Bantu language. From now on it
is also an open-source object-oriented library which has the ambition to be
generic and efficient. Even though the code is written to be generic, Akantu
strength are in solid mechanics models for fracture and contact simulations.

The full documentation can be found on [ReadTheDocs](https://akantu.readthedocs.io/en/latest/)

# Building `Akantu`

## Dependencies

In order to compile `Akantu`  any compiler supporting fully C++14 should work.
In addition some libraries are required:

 - CMake (>= 3.5.1)
 - Boost (pre-processor and Spirit)
 - Eigen3 (if not present the build system will try to download it)

For the python interface:

 - Python (>=3 is recommended)
 - pybind11 (if not present the build system will try to download it)

To run parallel simulations:

 - MPI
 - Scotch

To use the static or implicit dynamic solvers at least one of the following libraries is needed:

 - MUMPS (since this is usually compiled in static you also need MUMPS dependencies)
 - PETSc

To compile the tests and examples:

 - Gmsh
 - google-test (if not present the build system will try to download it)

### On `.deb` based systems

``` sh
 > sudo apt install cmake libboost-dev gmsh libeigen3-dev
 # For parallel
 > sudo apt install mpi-default-dev libmumps-dev libscotch-dev
 # For sequential
 > sudo apt install libmumps-seq-dev 
```

### Using `conda`

This works only for sequential computation since `mumps` from conda-forge is compiled without MPI support

``` sh
 > conda create -n akantu
 > conda activate akantu
 > conda install boost cmake
 > conda install -c conda-forge mumps
```

### Using `homebrew`

``` sh
 > brew install gcc boost scotch eigen
 > brew tap brewsci/num
 > brew install brewsci-mumps --without-brewsci-parmetis
```

If it does not work you can point "url" to "http://graal.ens-lyon.fr/MUMPS/MUMPS_5.6.2.tar.gz" together with "sha256" to "13a2c1aff2bd1aa92fe84b7b35d88f43434019963ca09ef7e8c90821a8f1d59a" using the command::

``` sh
  > brew edit brewsci/mumps
```

For the brew installation of mumps to work I also had to change in brewsci/num "alllib" -> "all", "-shared" -> "-dynamiclib" and remove all mentions of "-noall_load".

## Configuring and compilation


`Akantu` is a [CMake](https://cmake.org/) project, so to configure it, you can follow the usual way:

``` sh
  > cd akantu
  > mkdir build
  > cd build
  > ccmake ..
  [ Set the options that you need ]
  > make
  > make install

```

### On Mac OS X with `homebrew`

You will need to specify the compiler explicitly:

``` sh
 > CC=gcc-12 CXX=g++-12 FC=gfortran-12 cmake ..
```

Considering the homebrew is installed in `/opt/homebrew`
Define the location of the ``Scotch`` library path:
``` sh
 > cmake .. -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"
```

Specify path to all ``MUMPS`` libraries:
``` sh
 > cmake .. -DMUMPS_DIR=/opt/homebrew/opt/brewsci-mumps
```

In case the above does not work, specify the ``MUMPS`` path manually using (e.g.):
``` sh
 > cmake .. -DMUMPS_LIBRARY_COMMON=/opt/homebrew/opt/brewsci-mumps/lib/libmumps_common.dylib 
```

If compilation does not work change the path of the failing libraries to brew downloads in `/opt/homebrew/`. 

## Using the python interface


You can install ``Akantu`` using pip, this will install a pre-compiled version, this works only on Linux machines for now:

``` sh
  > pip install akantu
```

You can then import the package in a python script as:

``` python
  import akantu
```

The python API is similar to the C++ one. If you
encounter any problem with the python interface, you are welcome to do a merge
request or post an issue on [GitLab](https://gitlab.com/akantu/akantu/-/issues).

# Contributing

## Contributing new features, bug fixes

Any contribution is welcome, we are trying to follow a [gitflow](https://nvie.com/posts/a-successful-git-branching-model/) workflow, so the project `developers` can create branches named `features/<name of my feature>` or `bugfixes/<name of the fix>` directly in the main `akantu` repository.
External fellows can [Fork](https://gitlab.com/akantu/akantu/-/forks/new) the project.
In both cases the modifications have to be submitted in the form of a [Merge Request](https://gitlab.com/akantu/akantu/-/merge_requests/new).

## Asking for help, reporting issues

If you want to ask for help concerning Akantu's compilation, usage or problem with the code do not hesitate to open an [Issue](https://gitlab.com/akantu/akantu/-/issues/new) on gitlab. If you want to contribute and don't know where to start, you are also invited to open an issue.


# Examples and Tutorials with the python interface

To help getting started, you can find examples with the source code in the
`examples` sub-folder. If you just want to test the python examples without
having to compile the whole project you can use the following tarball
[akantu-python-examples.tgz](https://gitlab.com/akantu/akantu/-/packages/22034181).

In addition to the examples, multiple tutorials using the python interface are
available as notebooks with pre-installed version of `Akantu` on Renku. The
tutorials can be tested here:

[![renku](https://user-content.gitlab-static.net/52a4794df1236b248c8fc870bd74e9d787c0e2cb/68747470733a2f2f72656e6b756c61622e696f2f72656e6b752d62616467652e737667)](https://renkulab.io/projects/guillaume.anciaux/akantu-tutorials/sessions/new?autostart=1)

