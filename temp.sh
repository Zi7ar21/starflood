#!/bin/bash

# stop script if a command is unsuccessful
set -e

RED='\033[0;31m'
ORANGE="\033[0;33m"
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m'

# 0 - Normal Style
# 1 - Bold
# 2 - Dim
# 3 - Italic
# 4 - Underlined
# 5 - Blinking
# 7 - Reverse
# 8 - Invisible

mkdir -p out
echo -e "${RED}This script will mount a ${ORANGE}2 GiB${RED} tmpfs at \"./out\".${NC}" && \
echo -e "${ORANGE}tmpfs-es are meant to be fast ephemeral (non-volatile) storage, meaning the files inside will dissapear when you unmount it OR restart your machine.${NC}" && \
echo -e "When you are done, ${GREEN}you can unmount the tmpfs with \"${CYAN}sudo umount out${GREEN}.\"${NC}" && \
echo -e "" && \
echo -e "Press [${CYAN}ENTER${NC}] if you wish to continue" && \
echo -e "[${CYAN}CTRL${NC}]+[${CYAN}C${NC}] if you wish to signal SIGINT (exit)" && \
read -p ""

echo -e "Mounting a 2 GiB tmpfs at ./out..."

sudo mount -t tmpfs -o size=2g tmpfs out

echo -e "tmpfs mounted. Remember, ${GREEN}to unmount the tmpfs when you are done, you can use \"${CYAN}sudo umount out${GREEN}\".${NC}"
