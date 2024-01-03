#!/bin/bash

# stop script if a command is unsuccessful
set -e

echo "+-----------------------+" && \
echo "| Starflood Test Script |" && \
echo "+-----------------------+" && \
echo "" && \
echo "Cleaning build..." && \
echo "" && \
# todo
echo "Compiling..." && \
echo "" && \
./build.sh && \
echo "" && \
echo "Cleaning render..." && \
echo "" && \
rm -rf ./out/* && \
mkdir -p ./out && \
echo "Running..." && \
echo "" && \
/usr/bin/time -v build/starflood.out && \
echo "" && \
echo "Encoding render.mp4..." && \
rm -f ./render.mp4 && \
./encode.sh && \
echo "" && \
echo "Finished!"
