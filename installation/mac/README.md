# macOS Installation Package

The macOS installation package provides three command files that serve different purposes:

1. An installation file `install_fortran.command` to install Fortran on your Mac.
2. A file `uninstall_fortran.command` that completely removes the Fortran installation from your Mac.
3. A file `update_toolbox.command` that compiles and copie the latest version of our toolbox to your system.

All installation files have been tested with macOS Big Sur, but should work on previous versions.

## Important security motice for macOS users ##

Starting with macOS Catalina, Apple has removed terminal access to the "Documents" and the "Download" folder for security purposes. As a result, when you download an installation files from this respository (into Documents or Downloads) and try to execute the respective install_fortran.command file, it will not execute correctly. To solve this issue, please copy the installation files into another folder of yours, e.g. create a folder "Programming" in your home directory. If you install fortran from there, installation (and uninstallation) should run smoothly.

The same is true when compiling, building and executing Fortran (.f90) files. You can't do this within your Documents or Download folder anymore. Please copy your source files to another folder, again e.g. create a folder "Programming" in your home directory.

There are some fixes for this issue on the internet. Yet, they involve changing access privileges to your hard drive, which we can not recommend. If anyone still faces problems with this, please open a thread on the forum.

## Installation

When you execute the file `install_fortran.command`, it will install the [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) environment, the plotting programm [gnuplot](http://gnuplot.info/) as well as the IDE [geany](https://www.geany.org/) to your Mac. In addition, it is going to compile and store the toolbox on your Mac.


#### The installation process in detail

The file `install_fortran.command` is a so-called bash-file, which we use to execute certain installation routines on your Mac. You can run this file by right-clicking on it and then selecting `Open`. During the installation process, we will make some changes to your lower system files. For this, the file `install_fortran.command` needs administrator rights, which is why you are asked for your password twice. When starting the file, you will be asked to confirm that you want to install all software components as well as our toolbox on your Mac. You can do this by typing the letter `y` and pressing ENTER.

Now the installation process has started. In detail you will see the following things happening:

1. The first part of the installer copies all [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) components to your system using [homebrew](https://brew.sh).
2. Next, we install [gnuplot](http://gnuplot.info/).
3. The third installer component will bring the IDE [geany](https://www.geany.org/) to your system.

When the three software components are installed, the installation bash-file will do a couple more things:

1. It will patch the software geany such that it already contains the right configuration for using Fortran with our toolbox.
2. It will compile and install the toolbox to a directory to which the compiler is pointed.

**This is it!** If you saw no error messages appearing, then Fortran as well as our toolbox were successfully installed on your system.


## Uninstallation

In the unlikely event that you want to uninstall all the software components as well as the toolbox from your system, you can use the file `uninstall_fortran.command`. It follows the same logic as the installation bash file. You just have to follow the instructions.


## Toolbox Update

Sometimes you might feel that your toolbox is out of date. In this case, just pull the most recent version of this github repository and execute the file `update_toolbox.command`. This will compile the latest toolbox version and store is such that you can use it in all your programs.


## Run your first Fortran program

If you've successfully installed Fortran on you , [you can learn here how to run your first program](https://www.ce-fortran.com/run-program-mac/).
