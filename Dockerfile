FROM ubuntu:latest AS builder

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

# Install Movi
RUN git clone https://github.com/mohsenzakeri/Movi /movi && \
    cd /movi && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j

# Add the build directories to the PATH
ENV PATH="/movi/build:${PATH}"
