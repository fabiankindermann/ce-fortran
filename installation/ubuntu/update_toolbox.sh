#!/bin/bash

# A SHELL SCRIPT FOR UPDATING THE TOOLBOX ON UBUNTU
#
# ATTENTION: Fortran must already be installed using our original installation files.
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)


# set the current directory as running directory
cd "$( cd "$( dirname "$0" )" && pwd )"


# CHECK WHETHER THE SCRIPT HAS ROOT PRIVILIGES
[ "$UID" -eq 0 ] || { echo ; echo "THIS SCRIPT NEEDS TO BE RUN WITH ROOT PRIVILEGES!!!" ; echo ; echo "If you don't want this, use our docker images." ; echo ; echo "PLEASE TYPE YOUR PASSWORD:"; exec sudo "$0" "$@";}


# ASK FOR INSTALLATION CONFIRMATION
echo
echo "This script installs Fortran to your system."
echo
echo "ATTENTION: Fortran must already be installed using our original installation files."
echo
echo "THIS SCRIPT NEEDS ROOT PRIVILEGES FOR SOME INSTALLATION STEPS!!!"
echo 
read -rsp $'Do you want to continue (y/n)?' -n 1 key
echo

if [ "$key" != "y" ]; then
    exit 0
fi


## INSTALL THE TOOLBOX

# compile the toolbox
gfortran -c -Werror -fopenmp -Wno-unused -ffree-line-length-none -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -frecursive -g ./../toolbox/toolbox.f90 -o toolbox_debug.o
gfortran -c -O3 -fopenmp -ffree-line-length-none ./../toolbox/toolbox.f90 -o toolbox.o

# copy the toolbox to the working directory
sudo mkdir -p /usr/local/include
sudo mv toolbox.mod /usr/local/include/
sudo mv toolbox.o /usr/local/include/
sudo v toolbox_debug.o /usr/local/include/


# IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE INSTALLED PROPERLY
echo  
echo ...TOOLBOX UPDATE COMPLETED.
echo 
echo
echo In case you encountered any problem, check on www.ce-fortran.com for help.
echo
echo
read -p "Press RETURN to end..."
