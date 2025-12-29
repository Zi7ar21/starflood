#!/bin/sh

OUTPUT_DIR="./out"

FRAMERATE=30



#OPTIONS=-loop 0 -vf "scale=in_transfer=linear:out_transfer=bt709" -framerate ${FRAMERATE}
OPTIONS=-loop 0 -framerate ${FRAMERATE}

#INPUT_FILE="${OUTPUT_DIR}/vis/step_%04d.pfm"
#INPUT_FILE="${OUTPUT_DIR}/vis/step_%04d.ppm"
INPUT_FILE="${OUTPUT_DIR}/vis/step_%04d.png"

ffplay ${OPTIONS} -i ${INPUT_FILE}
