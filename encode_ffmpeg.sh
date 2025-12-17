#!/bin/sh

DATETIME=$(date +%Y%m%d_%H%M%S%N)

echo ${DATETIME}

ffmpeg -y -r 25 -framerate 25 -i ./out/vis/step_%04d.pfm -c:v h264_nvenc -vf "scale=in_transfer=linear:out_transfer=bt709" -b:v 8M -pix_fmt yuv420p -an -sn -dn -g 15 -bf 2 -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle starflood_${DATETIME}.mp4 && \
mpv --loop starflood_${DATETIME}.mp4
