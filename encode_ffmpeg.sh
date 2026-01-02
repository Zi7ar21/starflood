#!/bin/bash

# exit on error
set -e

OUTPUT_DIR="./out"
#OUTPUT_DIR="/scratch/Simulations/starflood_expanding_sphere_0"

FRAMERATE="30"



OPTIONS="-y"

INFILE_OPTIONS="-r ${FRAMERATE} -framerate ${FRAMERATE}"

#INFILE="${OUTPUT_DIR}/vis/step_%04d.pfm"
#INFILE="${OUTPUT_DIR}/vis/step_%04d.ppm"
INFILE="${OUTPUT_DIR}/vis/step_%04d.png"

#OUTFILE_OPTIONS="-c:v libx264 -b:v 8M -pix_fmt yuv420p -an -sn -dn -preset slow"
#OUTFILE_OPTIONS="-c:v h264_mediacodec -b:v 8M -pix_fmt nv12 -an -sn -dn"
#OUTFILE_OPTIONS="-c:v h264_nvenc -vf \"scale=in_transfer=linear:out_transfer=bt709\" -b:v 8M -pix_fmt yuv420p -an -sn -dn -g 15 -bf 2 -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle"
OUTFILE_OPTIONS="-c:v h264_nvenc -b:v 8M -pix_fmt yuv420p -an -sn -dn -g 15 -bf 2 -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle"

#OUTFILE="${OUTPUT_DIR}/starflood_vis_$(date +%Y%m%d_%H%M%S%N).mp4"
OUTFILE="${OUTPUT_DIR}/starflood_vis_$(date +%Y%m%d).mp4"
#OUTFILE="${OUTPUT_DIR}/starflood_vis.mp4"

echo "Encoding ${OUTFILE}..."

# Encode the visualization as a video file
ffmpeg ${OPTIONS} ${INFILE_OPTIONS} -i ${INFILE} ${OUTFILE_OPTIONS} ${OUTFILE}

echo "Finished encoding ${OUTFILE}."
