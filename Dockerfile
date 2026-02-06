FROM ubuntu:latest AS builder
ENV DEBIAN_FRONTEND=noninteractive
ENV CC=gcc-10
ENV CXX=g++-10
ENV CFLAGS="-march=x86-64 -mtune=generic -O2"
ENV CXXFLAGS="-march=x86-64 -mtune=generic -O2"

ARG MOVI_VERSION=v2.0.0

RUN apt-get update && apt-get install -y \
    zlib1g-dev \
    git \
    cmake \
    build-essential \
    python3 \
    wget \
    gcc-10 \
    g++-10 \
    time && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 10 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 10 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN git clone --branch ${MOVI_VERSION} --depth 1 https://github.com/mohsenzakeri/Movi /movi && \
    cd /movi && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
          -DCMAKE_C_FLAGS="${CFLAGS}" \
          -DCMAKE_BUILD_TYPE=Release \
          .. && \
    make -j$(nproc)

RUN cp /movi/build/bin/movi-* /usr/local/bin/

ENV PATH="/movi/build:/movi/build/bin:${PATH}"
