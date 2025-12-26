#!/bin/sh

make clean && \
make -j$(nproc) all && \
nice -n 1 build/starflood
