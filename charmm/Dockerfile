FROM ubuntu:24.04 AS build_charmm_base
RUN \
    apt-get update -y && \
    apt-get install -y libgfortran5 libgomp1

FROM build_charmm_base AS build_charmm
ADD charmm.tar.gz /

# Lines removing references to fstore.F90 are a workaround for a problem
# preventing CHARMM c51b1 from building.  The OpenMM interface for CHARMM in
# this release is broken.  This has nothing to do with OpenMMForceFields and
# this interface is not used in the force field conversion or testing process.
# Since OpenMM is not installed in the Docker image, the OpenMM interface will
# not be built so removing this file should have no effect, but its presence
# causes the build to fail for unknown reasons even if the OpenMM interface is
# explicitly disabled.
RUN \
    apt-get update -y && \
    apt-get install -y build-essential cmake gfortran && \
    cd /charmm && \
    rm source/openmm/fstore.F90 && \
    sed -i "/fstore/d" CMakeLists.txt && \
    sed -i "/fstore/d" tool/cmake/fortran_files.cmake && \
    ./configure --prefix /usr/local && \
    cd build/cmake && \
    make -j install

FROM build_charmm_base
COPY --from=build_charmm /usr/local/bin/charmm /usr/local/bin/charmm
COPY --from=build_charmm /usr/local/lib/libchmm.so /usr/local/lib/libchmm.so
