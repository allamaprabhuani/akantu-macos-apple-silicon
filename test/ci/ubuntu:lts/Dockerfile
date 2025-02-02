FROM ubuntu:22.04

ENV TZ=Europe/Zurich
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install akantu dependencies
RUN apt-get -qq update && apt-get -qq -y install \
    g++ gfortran  cmake gfortran g++\
    openmpi-bin libmumps-dev libscotch-dev \
    libboost-dev libopenblas-dev libcgal-dev \
    libeigen3-dev \
    python3 python3-dev python3-numpy python3-scipy python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install mpi4py

# Install test utilities
RUN apt-get -qq update && apt-get -qq -y install \
    python3-yaml python3-pytest python3-termcolor python3-click \
    ccache  gmsh curl git xsltproc \
    gcovr binutils \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install warning_parser

COPY .openmpi /root/.openmpi

# for debug
RUN apt-get -qq update && apt-get -qq -y install \
    gdb valgrind emacs-nox \
    && rm -rf /var/lib/apt/lists/*
