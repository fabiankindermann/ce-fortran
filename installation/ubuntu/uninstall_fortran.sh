#!/bin/bash

# A SHELL SCRIPT FOR UNINSTALLING FORTRAN FROM UBUNTU
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)


# set the current directory as running directory
cd "$( cd "$( dirname "$0" )" && pwd )"


# ASK FOR UNINSTALLATION CONFIRMATION
echo
echo "This script uninstalls Fortran/GNU Plot/Geany from your system."
echo
echo "THIS SCRIPT NEEDS ROOT PRIVILEGES FOR MANY UNINSTALLATION STEPS!!!"
echo "PLEASE USE WITH CAUTION!"
echo 
read -rsp $'Do you want to continue (y/n)?' -n 1 key
echo

if [ "$key" != "y" ]; then
    exit 0
fi


# REMOVE THE TOOLBOX
sudo rm -f /usr/local/include/toolbox.mod
sudo rm -f /usr/local/include/toolbox.o
sudo rm -f /usr/local/include/toolbox_debug.o
sudo rm -f /usr/local/include/toolbox_version.sh


# UNINSTALL GEANY

# uninstall software
sudo apt-get --yes remove geany
sudo rm -r ~/.config/geany

# remove desktop icon
sudo rm -f ~/Desktop/geany.desktop


# UNINSTALL GNUPLOT
sudo apt-get --yes remove gnuplot gnuplot-x11

# UNINSTALL GNU FORTRAN COMPILER
sudo apt-get --yes remove gfortran

# DELETE DEPENDENCIES
sudo apt-get --yes autoremove

# CLEAN UP CONFIGURATIONS
sudo sudo apt-get --yes clean


# IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE UNINSTALLED PROPERLY
echo  
echo ...UNINSTALLATION PROCESS COMPLETED.
echo 
echo
echo In case you encountered any problem, check on www.ce-fortran.com for help.
echo
echo
read -p "Press RETURN to end..."
