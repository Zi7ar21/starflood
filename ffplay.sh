#!/bin/sh

#ffplay -loop 0 -vf "scale=in_transfer=linear:out_transfer=bt709" -framerate 30 -i ./out/vis/step_%04d.pfm
ffplay -loop 0 -framerate 30 -i ./out/vis/step_%04d.png
