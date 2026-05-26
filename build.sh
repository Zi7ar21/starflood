#!/usr/bin/env sh

make clean && mkdir -p ./build && make -j$(nproc) all
