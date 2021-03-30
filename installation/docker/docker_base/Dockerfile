# A DOCKER FILE TO BUNDLE THE CE-FORTRAN ENVIRONMENT INTO A DOCKER CONTAINER WITH GUI
#
# Special thanks to Lars Vilhuber for starting this project.
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Fabian Kindermann (contact@ce-fortran.com)

FROM ubuntu:20.04

WORKDIR /home/user

# COPY SOURCEy FILES INTO BASIC DIRECTORY
COPY install_fortran.sh .

# COPY CE-FORTRAN PROGRAMS
COPY ./code-book ./code-book

# COPY CE-FORTRAN SOLUTIONS
COPY ./code-solution ./code-solution

# COPY TOOLBOX
COPY toolbox.f90 .

# ADD FONT FILE
COPY ./localfonts.conf /etc/fonts/local.conf

# INSTALL FORTRAN DEVELOPMENT ENVIRONMENT
RUN chmod a+rx install_fortran.sh \  
    && DEBIAN_FRONTEND="noninteractive" TZ="Europe/Berlin" ./install_fortran.sh \
    && rm -r install_fortran.sh && rm -rf /var/lib/apt/lists/*

# SET ENVIRONMENT VARIABLE FOR COMPILATIONS
ENV FCOMPILE="-O3 -fopenmp -frecursive -ffree-line-length-none -Wno-unused -fimplicit-none -Wall -fcheck=bound,do -ffpe-trap=invalid,zero,overflow -I/usr/local/include/ /usr/local/include/toolbox.o"


