#!/bin/bash

# A SHELL SCRIPT FOR INSTALLING FORTRAN TO A UBUNTU LINUX SYSTEM
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)


# INSTALL BUILD ESSENTIAL TOOLS IF NOT YET DONE
apt update
apt upgrade
apt --yes install build-essential


# INSTALL MS CORE FONTS

# this part is adapted from Frank Hoffsümmer (https://github.com/captnswing/msttcorefonts)
apt update
apt install -y --no-install-recommends software-properties-common curl
apt-add-repository multiverse
apt update

# ms core fonts
echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections
apt install -y --no-install-recommends fontconfig ttf-mscorefonts-installer
RUN fc-cache -f -v


# GFORTRAN COMPILER
apt --yes install gfortran

# GDB FOR DEBUGGING
apt --yes install gdb

# GNUPLOT
apt --yes install gnuplot


## INSTALL THE TOOLBOX

# compile the toolbox
gfortran -c -Werror -Wno-unused -ffree-line-length-none -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -frecursive -g ./toolbox.f90 -o toolbox_debug.o
gfortran -c -O3 -ffree-line-length-none ./toolbox.f90 -o toolbox.o

# copy the toolbox to the working directory
mkdir -p /usr/local/include
mv toolbox.mod /usr/local/include/
mv toolbox.o /usr/local/include/
mv toolbox_debug.o /usr/local/include/
echo "no" >> /usr/local/include/toolbox_visual.txt