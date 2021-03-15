#!/bin/bash


# A SHELL SCRIPT FOR UNINSTALLING FORTRAN FROM MAC OS X
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
# #VC# VERSION: 1.5  (06 January 2021)


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

# and copy the toolbox to the preferred working directory
echo
echo ...REMOVE TOOLBOX...
sudo rm -f /usr/local/include/toolbox.mod
sudo rm -f /usr/local/include/toolbox.o
sudo rm -f /usr/local/include/toolbox_debug.o
sudo rm -f /usr/local/include/toolbox_version.command
echo ...DONE...
echo

# completely uninstall homebrew with all packages
echo 
echo ...UNINSTALLING HOMEBREW...
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall.sh)"
echo ...DONE...
echo

# remove geany application
echo
echo ...UNINSTALLING GEANY...
rm -rf /Applications/Geany.app/
rm -rf ~/.config/geany/
echo ...DONE...
echo


# if everything ran correctly, at this point everything should be uninstalled properly
echo  
echo ...UNINSTALLATION PROCESS COMPLETED.
echo 
echo
echo In case you encountered any problem, check on www.ce-fortran.com for help.
echo
echo
read -rsp $'Press RETURN to end...\n' -n 1 key