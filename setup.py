#!/usr/bin/env python
# -*- coding: utf-8 -*-
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import sys
import os.path
import pybind11 as py11
import configparser
from setuptools import find_packages
from packaging.version import Version
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version

try:
    from skbuild import setup
except ImportError:
    sys.stderr.write(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in"
        " pyproject.toml yourself"
    )
    raise

# This is needed for semver.py to be importable
source_folder = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(source_folder, "cmake"))

parser = configparser.ConfigParser()
parser.read("setup.cfg")
cmake_args = ["-Dpybind11_DIR:PATH={}".format(py11.get_cmake_dir())]

if "cmake_config" in parser:
    for k, v in parser["cmake_config"].items():
        k = k.upper()
        cmake_args.append("-D{}:BOOL={}".format(k, v))

akantu_libs = []
if "CI_AKANTU_INSTALL_PREFIX" in os.environ:
    ci_akantu_install_prefix = os.environ["CI_AKANTU_INSTALL_PREFIX"]
    akantu_dir = os.path.join(ci_akantu_install_prefix,
                              "lib", "cmake", "Akantu")
    akantu_libs.extend(
        [
            # paths comming from the manylinux install via gitlab-ci
            "/softs/view/lib/*",
            "/softs/view/lib64/*",
            os.path.join(ci_akantu_install_prefix, "lib64/*"),
            os.path.join(ci_akantu_install_prefix, "lib/*"),
        ]
    )
    cmake_args.extend(
        [
            "-DAKANTU_BYPASS_AKANTU_TARGET:BOOL=ON",
            "-DAkantu_DIR:PATH={}".format(akantu_dir),
        ]
    )

setup_kw = {}
try:
    import semver

    _version = semver.get_version()
    setup_kw = {
        "version": _version,
    }
    cmake_args.append("-DAKANTU_VERSION={}".format(_version))
except ImportError:
    pass


# Add CMake as a build requirement if cmake is not installed or is too low a
# version
setup_requires = []
try:
    if Version(get_cmake_version()) < Version("3.4"):
        setup_requires.append("cmake")
except SKBuildError:
    setup_requires.append("cmake")

with open(os.path.join(source_folder, "README.md"), "r") as fh:
    long_description = fh.read()

setup(
    name="akantu",
    url="https://akantu.ch",
    author="Nicolas Richart",
    author_email="nicolas.richart@epfl.ch",
    description="Akantu: Swiss-Made Open-Source Finite-Element Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    platforms="",
    license="L-GPLv3",
    license_files=["COPYING", "COPYING.lesser"],
    project_urls={
        "Bug Tracker": "https://gitlab.com/akantu/akantu/-/issues",
        "Source": "https://gitlab.com/akantu/akantu",
    },
    setup_requires=setup_requires,
    install_requires=["numpy", "scipy"],
    package_data={"AkantuLibs": akantu_libs},
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    include_package_data=False,
    cmake_args=cmake_args,
    cmake_languages=["CXX"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",  # noqa
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Topic :: Education",
        "Topic :: Scientific/Engineering",
    ],
    **setup_kw,
)
