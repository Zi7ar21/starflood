#!/bin/sh

make clean && mkdir -p ./build && make -j$(nproc) all
