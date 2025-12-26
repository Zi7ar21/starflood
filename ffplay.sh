#!/bin/sh

FRAMERATE=30

#ffplay -loop 0 -vf "scale=in_transfer=linear:out_transfer=bt709" -framerate ${FRAMERATE} -i ./out/vis/step_%04d.pfm
#ffplay -loop 0 -framerate ${FRAMERATE} -i ./out/vis/step_%04d.ppm
ffplay -loop 0 -framerate ${FRAMERATE} -i ./out/vis/step_%04d.png
