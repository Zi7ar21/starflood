#ffmpeg -r 30 -i out/frames/img_%03d.png -crf 17 -framerate 30 -vcodec libx264 -pix_fmt yuv420p -preset slow render.mp4
#ffmpeg -framerate 30 -stats -r 30 -an -sn -i out/frames/img_%03d.png -c h264 -r 30 -an -sn -pix_fmt yuv420p render.mp4
/usr/bin/time -v ffmpeg -framerate 15 -i out/IMG_%04d.png -c:v libx264 -preset veryslow -crf 17 -an -sn -pix_fmt yuv420p render.mp4
#ffmpeg -i input -c:v libx264 -preset slow -crf 22 -c:a copy output.mkv
