#!/bin/sh

#LOG_FILENAME="starflood_$(date +%Y%m%d_%H%M%S%N).log"
LOG_FILENAME="starflood_$(date +%Y%m%d).log"

make clean && \
make -j$(nproc) all && \
nice -n 1 build/starflood \
#| tee ${LOG_FILENAME}
