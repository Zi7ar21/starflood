echo "+-----------------------+" && \
echo "| Starflood Test Script |" && \
echo "+-----------------------+" && \
echo "" && \
echo "Cleaning build..." && \
echo "" && \
# todo
echo "Compiling..." && \
echo "" && \
g++ -fopenmp -march=native -Ofast -pedantic -std=c++11 -Wall -Wextra -Wshadow -o build/starflood.out src/rng.c src/image.c src/main.cpp && \
#/usr/bin/time -v g++ -fopenmp -march=x86-64 -O0 -pedantic -std=c++11 -Wall -Wextra -Wshadow -o build/starflood.out src/rng.c src/image.c src/main.cpp && \
echo "" && \
echo "Cleaning render..." && \
echo "" && \
rm -rf ./out/* && \
mkdir -p ./out && \
rm -f ./render.mp4 && \
echo "Running..." && \
echo "" && \
/usr/bin/time -v build/starflood.out && \
echo "" && \
echo "Encoding render.mp4..." && \
ffmpeg -framerate 10 -i out/img_%04d.hdr -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4 && \
#/usr/bin/time -v ffmpeg -framerate 10 -i out/img_%04d.hdr -c:v libx264 -preset slow -crf 17 -an -sn -pix_fmt yuv420p render.mp4 && \
echo "" && \
echo "Finished!"
