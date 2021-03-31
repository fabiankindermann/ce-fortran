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
apt --yes install dconf-editor dbus-x11
apt --yes install gnome-terminal


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


# INSTALL GNU GFORTRAN COMPILER
apt --yes install gfortran

# INSTALL GNU DEBUGGER
apt --yes install gdb

# INSTALL GNUPLOT
apt --yes install gnuplot gnuplot-x11

# INSTALL GEANY TEXT EDITOR
apt --yes install yaru-theme-icon
apt --yes install geany


# PATCH GEANY FOR USE WITH FORTRAN

# Fortran compilation file
echo '[build-menu]' >  "filetypes.fortran"
echo 'FT_00_LB=_Compile' >> "filetypes.fortran"
echo 'FT_00_CM=mkdir -p Build && gfortran -O3 -fopenmp -frecursive -ffree-line-length-none -Wno-unused -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -I/usr/local/include/ -J"./Build" -c "%f" -o "Build/%e.o"' >> "filetypes.fortran"
echo 'FT_00_WD=' >> "filetypes.fortran"
echo 'FT_01_LB=_Build' >> "filetypes.fortran"
echo 'FT_01_CM=gfortran -O3 -fopenmp -J"./Build" /usr/local/include/toolbox.o "Build/%e.o" -o Build/prog' >> "filetypes.fortran"
echo 'FT_01_WD=' >> "filetypes.fortran"
echo 'FT_02_LB=_Debug' >> "filetypes.fortran"
echo 'FT_02_CM=mkdir -p Build && gfortran -fopenmp -frecursive -ffree-line-length-none -Wno-unused -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -I/usr/local/include/ -J"./Build" -g "%f" -o "Build/%e.o" "/usr/local/include/toolbox_debug.o"' >> "filetypes.fortran"
echo 'FT_02_WD=' >> "filetypes.fortran"
echo 'EX_00_LB=_Run Program' >> "filetypes.fortran"
echo 'EX_00_CM=clear && "Build/prog"' >> "filetypes.fortran"
echo 'EX_00_WD=' >> "filetypes.fortran"
echo 'EX_01_LB=_Run Debugger'>> "filetypes.fortran"
echo 'EX_01_CM=clear && gdb "Build/%e.o"' >> "filetypes.fortran"
echo 'EX_01_WD=' >> "filetypes.fortran"

# geany editor behavior
echo '[geany]' >  "geany.conf"
echo 'indent_type=0' >>  "geany.conf"
echo 'pref_main_load_session=false' >>  "geany.conf"
echo 'beep_on_errors=false' >>  "geany.conf"

# update geany configuration
mkdir -p /root/.config/geany/filedefs/
mv ./filetypes.fortran /root/.config/geany/filedefs/filetypes.fortran
mv ./geany.conf /root/.config/geany/geany.conf


## INSTALL THE TOOLBOX

# compile the toolbox
gfortran -c -Werror -Wno-unused -ffree-line-length-none -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -frecursive -g ./toolbox.f90 -o toolbox_debug.o
gfortran -c -O3 -ffree-line-length-none ./toolbox.f90 -o toolbox.o

# copy the toolbox to the working directory
mkdir -p /usr/local/include
mv toolbox.mod /usr/local/include/
mv toolbox.o /usr/local/include/
mv toolbox_debug.o /usr/local/include/