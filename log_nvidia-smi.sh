#!/bin/sh

#FILENAME="./out/log_nvidia-smi_$(date +%Y-%m-%d).csv"
FILENAME="./out/log_nvidia-smi.csv"

#nvidia-smi --help-query-gpu
nvidia-smi --query-gpu timestamp,fan.speed,pstate,utilization.gpu,utilization.memory,temperature.gpu,power.draw,clocks.current.graphics,clocks.current.memory --format=csv,nounits --loop-ms=500 | tee ${FILENAME}
