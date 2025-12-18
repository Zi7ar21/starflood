#!/bin/sh

#FILENAME="starflood_$(date +%Y%m%d_%H%M%S%N).mp4"
FILENAME="starflood_2025-12-17.mp4"
#FILENAME="starflood_out.mp4"

echo ${FILENAME}

#ffmpeg -y -r 30 -framerate 30 -i ./out/vis/step_%04d.ppm -c:v h264_nvenc -vf "scale=in_transfer=linear:out_transfer=bt709" -b:v 16M -pix_fmt yuv420p -an -sn -dn -g 15 -bf 2 -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle ${FILENAME} && \
ffmpeg -y -r 30 -framerate 30 -i ./out/vis/step_%04d.ppm -c:v h264_nvenc -b:v 16M -pix_fmt yuv420p -an -sn -dn -g 15 -bf 2 -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle ${FILENAME} && \
mpv --loop ${FILENAME}
