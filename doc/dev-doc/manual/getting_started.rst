Getting Started
===============

Contributing
------------

Contributing new features, bug fixes
````````````````````````````````````

Any contribution is welcome, we are trying to follow a `gitflow <https://nvie.com/posts/a-successful-git-branching-model/>`_ workflow, so the project `developers` can create branches named `features/<name of my feature>` or `bugfixes/<name of the fix>` directly in the main `akantu` repository.
External fellows can `Fork <https://gitlab.com/akantu/akantu/-/forks/new>`_ the project.
In both cases the modifications have to be submitted in the form of a `Merge Request <https://gitlab.com/akantu/akantu/-/merge_requests/new>`_.

Asking for help, reporting issues
`````````````````````````````````

If you want to ask for help concerning Akantu's compilation, usage or problem with the code do not hesitate to open an `Issue <https://gitlab.com/akantu/akantu/-/issues/new>`_ on gitlab. If you want to contribute and don't know where to start, you are also invited to open an issue.


Building ``Akantu``
--------------------

Dependencies
````````````

In order to compile ``Akantu``  any compiler supporting fully C++14 should work.
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

On ``.deb`` based Linux systems
"""""""""""""""""""""""""""""""

.. code-block:: bash

 > sudo apt install cmake libboost-dev gmsh libeigen3-dev
 # For parallel
 > sudo apt install mpi-default-dev libmumps-dev libscotch-dev
 # For sequential
 > sudo apt install libmumps-seq-dev

Using ``conda``
"""""""""""""""

This works only for sequential computation since `mumps` from conda-forge is compiled without MPI support:

.. code-block:: bash

 > conda create -n akantu
 > conda activate akantu
 > conda install boost cmake
 > conda install -c conda-forge mumps

Using ``homebrew``
""""""""""""""""""

.. code-block:: bash

 > brew install gcc boost scotch eigen
 > brew tap brewsci/num
 > brew install brewsci-mumps --without-brewsci-parmetis

If it does not work you can point "url" to "http://graal.ens-lyon.fr/MUMPS/MUMPS_5.6.2.tar.gz" together with "sha256" to "13a2c1aff2bd1aa92fe84b7b35d88f43434019963ca09ef7e8c90821a8f1d59a" using the command::

.. code-block:: bash

 > brew edit brewsci-mumps

For the brew installation of mumps to work I also had to change in brewsci/num "alllib" -> "all", "-shared" -> "-dynamiclib" and remove all mentions of "-noall_load".

Configuring and compilation
```````````````````````````

`Akantu` is a `CMake <https://cmake.org/>`_ project, so to configure it, you can follow the usual way:

 .. code-block:: bash

  > cd akantu
  > mkdir build
  > cd build
  > ccmake ..
  [ Set the options that you need ]
  > make
  > make install

On Mac OS X with ``homebrew``
"""""""""""""""""""""""""""""
You will need to specify the compiler explicitly

.. code-block:: bash

  > CC=gcc-12 CXX=g++-12 FC=gfortran-12 cmake ..

Considering that `homebrew` is installed in ``/opt/homebrew``
Define the location of the ``Scotch`` library path:

.. code-block:: bash

 > cmake .. -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"

Specify path to all ``MUMPS`` libraries:

.. code-block:: bash

 > cmake .. -DMUMPS_DIR=/opt/homebrew/opt/brewsci-mumps

In case the above does not work, specify the ``MUMPS`` path manually using (e.g.):

.. code-block:: bash

 > cmake .. -DMUMPS_LIBRARY_COMMON=/opt/homebrew/opt/brewsci-mumps/lib/libmumps_common.dylib

If compilation does not work change the path of the failing libraries to brew downloads in `/opt/homebrew/`.

Using the python interface
--------------------------

You can install ``Akantu`` using pip, this will install a pre-compiled version, this works only on Linux machines for now::

  > pip install akantu

You can then import the package in a python script as::

  import akantu

The python API is similar to the C++ one, see :ref:`reference` . If you encouter any problem with the python interface, you are welcome to do a merge request or post an issue on `GitLab <https://gitlab.com/akantu/akantu/-/issues>`_ .
  

Examples and Tutorials with the python interface
````````````````````````````````````````````````
To help getting started, you can find examples with the source code in the
`examples` sub-folder. If you just want to test the python examples without
having to compile the whole project you can use the following tarball
`akantu-python-examples.tgz
<https://gitlab.com/akantu/akantu/-/packages/22034181>`_.

In addition to the examples, multiple tutorials using the python interface are
available as notebooks with pre-installed version of `Akantu` on `Renku`. The
tutorials can be tested here: |renku|

.. |renku| image:: https://user-content.gitlab-static.net/52a4794df1236b248c8fc870bd74e9d787c0e2cb/68747470733a2f2f72656e6b756c61622e696f2f72656e6b752d62616467652e737667
   :target: https://renkulab.io/projects/guillaume.anciaux/akantu-tutorials/sessions/new?autostart=1



Writing a ``main`` function
---------------------------

``Akantu`` first needs to be initialized. The memory management included in the
core library handles the correct allocation and de-allocation of vectors,
structures and/or objects. Moreover, in parallel computations, the
initialization procedure performs the communication setup. This is achieved by
the function :cpp:func:`initialize <akantu::initialize>` that is used as
follows::

    #include "aka_common.hh"
    #include "..."

    using namespace akantu;

    int main(int argc, char *argv[]) {
      initialize("input_file.dat", argc, argv);

      // your code ...

    }

The :cpp:func:`initialize <akantu::initialize>` function takes the text input
file and the program parameters which can be parsed by ``Akantu`` in due form
(see sect:parser). Obviously it is necessary to include all files needed in
main. In this manual, all provided code implies the usage of ``akantu`` as
namespace.

Compiling your simulation
-------------------------

The easiest way to compile your simulation is to create a ``cmake`` project by
putting all your code in some directory of your choosing. Then, make sure that
you have ``cmake`` installed and create a ``CMakeLists.txt`` file. An example of
a minimal ``CMakeLists.txt`` file would look like this:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.12.0)
   project(my_simu)

   find_package(Akantu REQUIRED)

   add_akantu_simulation(my_simu my_simu.cc)

Then create a directory called ``build`` and inside it execute ``cmake
-DAkantu_DIR=<path_to_akantu> -DCMAKE_BUILD_TYPE=RelWithDebInfo ..``. If you
installed ``Akantu`` in a standard directory such as ``/usr/local`` (using
``make install``), you can omit the ``-DAkantu_DIR=<path_to_akantu>`` option.

Otherwise ``path_to_akantu`` is either the folder where you built ``Akantu`` if
you did not do a ``make install``, or if you installed ``Akantu`` in
``CMAKE_INSTALL_PREFIX`` it is ``<CMAKE_INSTALL_PREFIX>/share/cmake/Akantu``.

Once ``cmake`` managed to configure and generate a ``makefile`` you can just do
``make``.


.. _loading_mesh:

Creating and Loading a Mesh
---------------------------

In its current state, ``Akantu`` supports three types of meshes: Gmsh, Abaqus and
Diana. Once a :cpp:class:`akantu::Mesh` object is created with a given spatial
dimension, it can be filled by reading a mesh input file. The method
:cpp:func:`read <akantu::Mesh::read>` of the class :cpp:class:`Mesh
<akantu::Mesh>` infers the mesh type from the file extension. If a non-standard
file extension is used, the mesh type has to be specified.

.. code-block:: c++

    Int spatial_dimension = 2;
    Mesh mesh(spatial_dimension);

    // Reading Gmsh files
    mesh.read("my_gmsh_mesh.msh");
    mesh.read("my_gmsh_mesh", _miot_gmsh);

The Gmsh reader adds the geometrical and physical tags as mesh data. The
physical values are stored as a :cpp:type:`Int <akantu::Int>` data called
``tag_0``, if a string name is provided it is stored as a ``std::string`` data
named ``physical_names``. The geometrical tag is stored as a :cpp:type:`Int
<akantu::Int>` data named ``tag_1``.

Running parallel simulation
---------------------------

In order to run distributed memory simulation a few extra steps have to be taken.
The mesh as to be distributed

.. code-block:: c++

    const auto & comm = Communicator::getStaticCommunicator();
    if (comm.whoAmI() == 0) {  // MPI rank
      // Read the mesh
      mesh.read("square_2d.msh");
    }
    mesh.distribute();

All the communications and the distribution of the mesh and associated data will
be taken care automatically.

Currently the mesh decomposition is handled by the `Scotch
<https://gitlab.inria.fr/scotch/scotch>`_ library. Which means if needed you
could define different edge and vertex weights

.. code-block:: c++

     mesh.distribute(_edge_weight_function =
                         [](auto &&, auto &&) { return 1; },
                      _vertex_weight_function =
                         [](auto &&) { return 1; });

The `vertex` weights correspond to the computational cost of the elements, and
the `edge` weights relates to the cost of communications between 2 elements.

To run the simulation you will need to use a runner appropriate to your machine,
like `mpirun`, `srun`, `arun`, etc.

.. code-block:: sh

  $ mpirun -np 4 ./my_simulation
