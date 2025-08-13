FROM ubuntu:latest AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set portable compiler flags to prevent illegal instruction errors
ENV CC=gcc-10
ENV CXX=g++-10
ENV CFLAGS="-march=x86-64 -mtune=generic -O2"
ENV CXXFLAGS="-march=x86-64 -mtune=generic -O2"

# Install required dependencies in a single layer
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

# Install Movi with portable compilation flags
RUN git clone https://github.com/mohsenzakeri/Movi /movi && \
    cd /movi && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
          -DCMAKE_C_FLAGS="${CFLAGS}" \
          -DCMAKE_BUILD_TYPE=Release \
          .. && \
    make -j$(nproc)

# TODO: Quick fix for finding hte relavant binaries
RUN cp /movi/build/movi-* /

# Add the build directories to the PATH
ENV PATH="/movi/build:${PATH}"
