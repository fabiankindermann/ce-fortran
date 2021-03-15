#!/bin/bash


# A SHELL SCRIPT FOR UNINSTALLING FORTRAN FROM MACOS
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)


# set the current directory as running directory
cd "$( cd "$( dirname "$0" )" && pwd )"



# ASK FOR UNINSTALLATION CONFIRMATION

echo
echo This script uninstalls Fortran/GNU Plot/Geany from your system.
echo 
read -rsp $'Do you want to continue (y/n)?' -n 1 key
echo

if [ "$key" != "y" ]; then
    exit 0
fi



# REMOVE THE TOOLBOX FILES

sudo rm -f /usr/local/include/toolbox.mod
sudo rm -f /usr/local/include/toolbox.o
sudo rm -f /usr/local/include/toolbox_debug.o
sudo rm -f /usr/local/include/toolbox_version.command



# COMPLETELY UNINSTALL HOMEBREW

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall.sh)"


# REMOVE GEANY APPLICATION

rm -rf /Applications/Geany.app/
rm -rf ~/.config/geany/



# IF EVERYTHING RAN CORRECTLY, AT THIS POINT EVERYTHING SHOULD BE UNINSTALLED PROPERLY

echo  
echo ...UNINSTALLATION PROCESS COMPLETED.
echo 
echo
echo In case you encountered any problem, check on www.ce-fortran.com for help.
echo
echo
read -rsp $'Press RETURN to end...\n' -n 1 key