FROM debian:bookworm

# library dependencies
RUN apt -qq update && apt -qq -y install \
    g++ gfortran clang cmake \
    openmpi-bin libmumps-dev libscotch-dev \
    libboost-dev libopenblas-dev \
    libeigen3-dev \
    python3 python3-dev python3-numpy python3-scipy python3-mpi4py \
    && rm -rf /var/lib/apt/lists/*

# for documentation
RUN apt -qq update && apt -qq -y install \
    python3-sphinx \
    python3-sphinxcontrib.bibtex \
    python3-sphinx-rtd-theme \
    python3-breathe \
    python3-git python3-jinja2 \
    doxygen graphviz \
    && rm -rf /var/lib/apt/lists/*

# for ci
RUN apt -qq update && apt -qq -y install \
    gmsh python3-pytest python3-click python3-termcolor \
    python3-flake8 python3-pip python3-tqdm \
    ccache clang-format clang-tidy \
    curl git xsltproc jq \
    gcovr llvm binutils \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages warning_parser

COPY .openmpi /root/.openmpi

# for debug
RUN apt -qq update && apt -qq -y install \
    gdb valgrind \
    && rm -rf /var/lib/apt/lists/*
