FROM quay.io/pypa/manylinux_2_28_x86_64 as builder

ENV LANG en_US.UTF-8

RUN yum install -y xz openssl-devel perl-devel
# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /softs/spack-environment
COPY spack.yaml /softs/spack-environment

RUN adduser -u 1000  -s /bin/false app
RUN mkdir -p /softs && chown -R app:app /softs

USER app

RUN cd /softs && git clone https://github.com/spack/spack.git -b v0.22.2

#RUN cd /softs && git clone https://gitlab.com/epfl-lsms/spack-packages.git

COPY packages.yaml /softs/spack/etc/spack
#COPY repos.yaml /softs/spack/etc/spack

# Install the software, remove unnecessary deps
RUN python3 /softs/spack/bin/spack -e /softs/spack-environment install && \
    python3 /softs/spack/bin/spack gc -y

# Strip all the binaries
RUN find -L /softs/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip -s

USER root

# Modifications to the environment that are necessary to run
RUN cd /softs/spack-environment && \
    python3 /softs/spack/bin/spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Bare OS image to run the installed executables
FROM quay.io/pypa/manylinux_2_28_x86_64

RUN yum install -y xz zlib-devel

RUN adduser -u 1000  -s /bin/false app

COPY --from=builder /softs/spack /softs/spack
COPY --from=builder /softs/spack-environment /softs/spack-environment
COPY --from=builder /softs/software /softs/software
COPY --from=builder /softs/view /softs/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

RUN chown -R app:app /softs

RUN mkdir -p /build && chown app:app /build

ENV CROSS_ROOT /opt/rh/gcc-toolset-13/root/usr/bin
ENV AS=${CROSS_ROOT}/as \
    AR=${CROSS_ROOT}/ar \
    CC=${CROSS_ROOT}/gcc \
    CPP=${CROSS_ROOT}/cpp \
    CXX=${CROSS_ROOT}/g++ \
    LD=${CROSS_ROOT}/ld \
    FC=${CROSS_ROOT}/gfortran \
    CMAKE_PREFIX_PATH=/soft/view

COPY install-pip.sh /build
RUN /build/install-pip.sh

COPY Toolchain.cmake ${CROSS_ROOT}/../lib/
ENV CMAKE_TOOLCHAIN_FILE ${CROSS_ROOT}/../lib/Toolchain.cmake

USER app

#ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]
