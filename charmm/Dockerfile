FROM ubuntu:24.04 AS build_charmm_base
RUN \
    apt-get update -y && \
    apt-get install -y libgfortran5 libgomp1

FROM build_charmm_base AS build_charmm
ADD charmm.tar.gz /

RUN \
    apt-get update -y && \
    apt-get install -y build-essential cmake gfortran ninja-build && \
    cd /charmm && \
    ./configure --with-ninja --prefix /usr/local && \
    cd build/cmake && \
    ninja install

FROM build_charmm_base
COPY --from=build_charmm /usr/local/bin/charmm /usr/local/bin/charmm
COPY --from=build_charmm /usr/local/lib/libchmm.so /usr/local/lib/libchmm.so
