#!/bin/bash


# A SHELL SCRIPT FOR INSTALLING FORTRAN ON UBUNTU
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)


# set the current directory as running directory
cd "$( cd "$( dirname "$0" )" && pwd )"



# ASK FOR INSTALLATION CONFIRMATION

echo
echo "This script installs Fortran to your system."
echo
echo "THIS SCRIPT NEEDS ROOT PRIVILEGES FOR MANY INSTALLATION STEPS!!!"
echo 
echo "If you don't want this, use our docker images." 
echo 
read -rsp $'Do you want to continue (y/n)?' -n 1 key
echo

if [ "$key" != "y" ]; then
    exit 0
fi



# INSTALL BUILD ESSENTIAL TOOLS IF NOT YET DONE

sudo apt-get update
sudo apt-get --yes install build-essential



# INSTALL GNU GFORTRAN COMPILER

sudo apt-get --yes install gfortran



# INSTALL GNU DEBUGGER

sudo apt-get --yes install gdb



# INSTALL GNUPLOT

sudo apt-get --yes install gnuplot gnuplot-x11



# INSTALL GEANY TEXT EDITOR

sudo apt-get --yes install geany



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
mkdir -p ~/.config/geany/filedefs/
cp ./src/filetypes.fortran ~/.config/geany/filedefs/filetypes.fortran
cp ./src/geany.conf ~/.config/geany/geany.conf
sudo chmod 664 ~/.config/geany/filedefs/filetypes.fortran
sudo chown $USER:$USER -c ~/.config/geany/filedefs/filetypes.fortran
sudo chmod 664 ~/.config/geany/geany.conf
sudo chown $USER:$USER -c ~/.config/geany/geany.conf

# geany desktop icon
cp /usr/share/applications/geany.desktop ~/Desktop
sudo chmod +x ~/Desktop/geany.desktop
sudo chown $USER:$USER -c ~/Desktop/geany.desktop


## INSTALL THE TOOLBOX

# compile the toolbox
gfortran -c -Werror -Wno-unused -ffree-line-length-none -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -frecursive -g ./../toolbox/toolbox.f90 -o toolbox_debug.o
gfortran -c -O3 -ffree-line-length-none ./../toolbox/toolbox.f90 -o toolbox.o

# copy the toolbox to the working directory
sudo sudo mkdir -p /usr/local/include
sudo mv toolbox.mod /usr/local/include/
sudo mv toolbox.o /usr/local/include/
sudo mv toolbox_debug.o /usr/local/include/



# IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE INSTALLED PROPERLY

echo  
echo ...INSTALLATION COMPLETED.
echo 
echo
echo In case you encountered any problem, check on www.ce-fortran.com for help.
echo
echo
read -p "Press RETURN to end..."
