FROM ubuntu:latest AS build_charmm_base
RUN \
    apt-get update -y && \
    apt-get install -y libgfortran5 libgomp1

FROM build_charmm_base AS build_charmm
ADD charmm.tar.gz /
RUN \
    apt-get update -y && \
    apt-get install -y build-essential cmake gfortran && \
    cd /charmm && \
    mkdir cmake && \
    cd cmake && \
    cmake .. && \
    make install

FROM build_charmm_base
COPY --from=build_charmm /usr/local/bin/charmm /usr/local/bin/charmm
