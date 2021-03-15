#!/bin/bash


# A SHELL SCRIPT FOR UNINSTALLING FORTRAN FROM UBUNTU
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)


# set the current directory as running directory
cd "$( cd "$( dirname "$0" )" && pwd )"



# CHECK WHETHER THE SCRIPT HAS ROOT PRIVILIGES

[ "$UID" -eq 0 ] || { echo ; echo "THIS SCRIPT NEEDS TO BE RUN WITH ROOT PRIVILEGES!!!" ; echo ; echo "PLEASE TYPE YOUR PASSWORD:"; exec sudo "$0" "$@";}



# ASK FOR UNINSTALLATION CONFIRMATION
echo
echo This script uninstalls Fortran/GNU Plot/Geany from your system.
echo 
read -rsp $'Do you want to continue (y/n)?' -n 1 key
echo

if [ "$key" != "y" ]; then
    exit 0
fi



# REMOVE THE TOOLBOX

rm -f /usr/local/include/toolbox.mod
rm -f /usr/local/include/toolbox.o
rm -f /usr/local/include/toolbox_debug.o
rm -f /usr/local/include/toolbox_version.sh


# UNINSTALL GEANY

# uninstall software
apt-get --yes remove geany
rm -r ~/.config/geany

# remove desktop icon
sudo rm -f ~/Desktop/geany.desktop



# UNINSTALL GNUPLOT

apt-get --yes remove gnuplot gnuplot-x11



# UNINSTALL GNU FORTRAN COMPILER

apt-get --yes remove gfortran



# DELETE DEPENDENCIES

apt-get --yes autoremove



# CLEAN UP CONFIGURATIONS

sudo apt-get --yes clean




# IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE UNINSTALLED PROPERLY

echo  
echo ...UNINSTALLATION PROCESS COMPLETED.
echo 
echo
echo In case you encountered any problem, check on www.ce-fortran.com for help.
echo
echo
read -p "Press RETURN to end..."
