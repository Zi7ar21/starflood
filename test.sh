echo "+-----------------------+" && \
echo "| Starflood Test Script |" && \
echo "+-----------------------+" && \
echo "" && \
echo "Cleaning..." && \
echo "" && \
rm -rf ./out && \
mkdir -p ./out && \
rm -f ./render.mp4 && \
echo "Compiling..." && \
echo "" && \
/usr/bin/g++ -O0 -g -pedantic -std=c++11 -Wall -Wextra -Wshadow -o build/starflood.out src/image.cpp src/main.cpp  && \
#/usr/bin/time -v g++ -Ieigen-3.4.0 -march=native -O2 -pedantic -std=c++11 -Wall -Wextra -Wshadow -o build/starflood.out src/image.cpp src/main.cpp && \
echo "" && \
echo "Running..." && \
echo "" && \
/usr/bin/time -v build/starflood.out && \
echo "" && \
echo "Encoding render.mp4..." && \
ffmpeg -framerate 15 -i out/IMG_%04d.tga -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4 && \
#/usr/bin/time -v ffmpeg -framerate 15 -i out/IMG_%04d.hdr -c:v libx264 -preset veryslow -crf 17 -an -sn -pix_fmt yuv420p render.mp4 && \
echo "" && \
echo "Finished!"
