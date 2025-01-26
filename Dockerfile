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

# Install pfp-thresholds
RUN git clone https://github.com/mohsenzakeri/pfp-thresholds  /pfp-thresholds && \
    cd /pfp-thresholds && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# Install r-permute
RUN git clone https://github.com/drnatebrown/r-permute /r-permute && \
    cd /r-permute && \
    git submodule update --init --recursive && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# Install Movi
RUN git clone https://github.com/mohsenzakeri/Movi /movi && \
    cd /movi && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# Configure preprocess_ref.sh
RUN sed -i 's|pfp=.*|pfp=/pfp-thresholds/build/pfp_thresholds|' /movi/preprocess_ref.sh && \
    sed -i 's|movi_default=.*|movi_default=/movi/build/movi-default|' /movi/preprocess_ref.sh && \
    sed -i 's|movi_constant=.*|movi_constant=/movi/build/movi-constant|' /movi/preprocess_ref.sh && \
    sed -i 's|prepare_ref=.*|prepare_ref=/movi/build/prepare_ref|' /movi/preprocess_ref.sh && \
    sed -i 's|bconstructor=.*|bconstructor=/r-permute/build/test/src/build_constructor|' /movi/preprocess_ref.sh && \
    sed -i 's|rconstructor=.*|rconstructor=/r-permute/build/test/src/run_constructor|' /movi/preprocess_ref.sh

# Add the build directories to the PATH
ENV PATH="/r-permute/build/test/src/:/pfp-thresholds/build:/movi/build:${PATH}"