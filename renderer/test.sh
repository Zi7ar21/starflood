echo "+-----------------------+" && \
echo "| Starflood Test Script |" && \
echo "+-----------------------+" && \
echo "" && \
echo "Cleaning build..." && \
echo "" && \
# todo
echo "Compiling..." && \
echo "" && \
g++ -fopenmp -march=native -Og -pedantic -std=c++17 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -o build/starflood.out src/main.cpp && \
#/usr/bin/time -v g++ -fopenmp -march=x86-64 -O0 -pedantic -std=c++11 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -o build/starflood.out src/main.cpp && \
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
# DO NOT TOUCH THIS COMMAND
# I SPENT HOURS FIGURING OUT
# HOW THE HECK TO GET FFMPEG
# TO RENDER LINEAR COLORS
# IN A SOMEWHAT CORRECT WAY
#ffmpeg -framerate 30 -colorspace bt709 -color_primaries bt709 -color_trc linear -i out/step_%04d.hdr -vf zscale=t=linear:npl=100,format=gbrpf32le,zscale=p=bt709,tonemap=tonemap=clip:desat=0,zscale=t=bt709:m=bt709:r=tv,format=yuv420p -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4 && \
#/usr/bin/time -v ffmpeg -framerate 30 -i out/step_%04d.hdr -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4 && \
echo "" && \
echo "Finished!"
