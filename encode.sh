#!/bin/bash

#ffmpeg -r 30 -i out/frames/img_%03d.png -crf 17 -framerate 30 -vcodec libx264 -pix_fmt yuv420p -preset slow render.mp4
#ffmpeg -framerate 30 -stats -r 30 -an -sn -i out/frames/img_%03d.png -c h264 -r 30 -an -sn -pix_fmt yuv420p render.mp4
#/usr/bin/time -v ffmpeg -framerate 15 -i out/IMG_%04d.png -c:v libx264 -preset veryslow -crf 17 -an -sn -pix_fmt yuv420p render.mp4
#ffmpeg -i input -c:v libx264 -preset slow -crf 22 -c:a copy output.mkv

# DO NOT TOUCH THIS COMMAND
# I SPENT HOURS FIGURING OUT
# HOW THE HECK TO GET FFMPEG
# TO RENDER LINEAR COLORS
# IN A SOMEWHAT CORRECT WAY
#ffmpeg -framerate 30 -colorspace bt709 -color_primaries bt709 -color_trc linear -i out/step_%04d.hdr -vf zscale=t=linear:npl=100,format=gbrpf32le,zscale=p=bt709,tonemap=tonemap=clip:desat=0,zscale=t=bt709:m=bt709:r=tv,format=yuv420p -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4
ffmpeg -framerate 30 -colorspace bt709 -color_primaries bt709 -color_trc linear -i out/step_%04d.hdr -vf zscale=t=linear:npl=100,format=gbrpf32le,zscale=p=bt709,tonemap=tonemap=clip:desat=0,zscale=t=bt709:m=bt709:r=tv,format=yuv420p -c:v h264_nvenc -preset slow -tune hq -b:v 8M -an -sn -pix_fmt yuv420p render.mp4
#/usr/bin/time -v ffmpeg -framerate 30 -i out/step_%04d.hdr -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4
