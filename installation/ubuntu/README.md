# Ubuntu Installation Package

The Ubuntu installation package provides three command files that serve different purposes:

1. An installation file `install_fortran.sh` to install Fortran on your computer.
2. A file `uninstall_fortran.sh` that completely removes the Fortran installation from your computer.
3. A file `update_toolbox.sh` that compiles and copie the latest version of our toolbox to your system.

All installation files have been tested with Ubuntu 20.04, but should work on previous and later versions.


## Installation

When you execute the file `install_fortran.sh`, it will install the [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) environment, the plotting programm [gnuplot](http://gnuplot.info/) as well as the IDE [geany](https://www.geany.org/) to your computer. In addition, it is going to compile and store the toolbox on your computer.


#### The installation process in detail

The file `install_fortran.sh` is a bash-file, which we use to execute certain installation routines on your computer. You can run this file in a terminal. You might have to add executable rights to it by typing `chmod +x install_fortran.sh`. During the installation process, we will make some changes to your lower system files. For this, the file `install_fortran.sh` needs administrator rights, which is why you are asked for your password. When starting the file, you will be asked to confirm that you want to install all software components as well as our toolbox on your computer. You can do this by typing the letter `y` and pressing ENTER.

Now the installation process has started. In detail you will see the following things happening:

1. The first part of the installer copies all [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) components to your system.
2. Next, we install [gnuplot](http://gnuplot.info/).
3. The third installer component will bring the IDE [geany](https://www.geany.org/) to your system.

When the three software components are installed, the installation bash-file will do a couple more things:

1. It will patch the software geany such that it already contains the right configuration for using Fortran with our toolbox.
2. It will compile and install the toolbox to a directory to which the compiler is pointed.

**This is it!** If you saw no error messages appearing, then Fortran as well as our toolbox were successfully installed on your system.


## Uninstallation

In the unlikely event that you want to uninstall all the software components as well as the toolbox from your system, you can use the file `uninstall_fortran.sh`. It follows the same logic as the installation bash file. You just have to follow the instructions.


## Toolbox Update

Sometimes you might feel that your toolbox is out of date. In this case, just pull the most recent version of this github repository and execute the file `update_toolbox.sh`. This will compile the latest toolbox version and store is such that you can use it in all your programs.


## Run your first Fortran program

If you've successfully installed Fortran on you , [you can learn here how to run your first program](https://www.ce-fortran.com/run-program-lin/).
