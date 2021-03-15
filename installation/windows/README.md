# Windows Installation Package

The Windows installation package provides three batch files that serve different purposes:

1. An installation file `install_fortran.bat` to install Fortran on your Windows computer.
2. A file `uninstall_fortran.bat` that completely removes the Fortran installation from your computer.
3. A file `update_toolbox.bat` that compiles and copies the latest version of our toolbox to your system.

All installation files have been tested with Windows 10, but should work down to Windows 7, as long as Power Shell is available.



## Installation

When you execute the file `install_fortran.bat`, it will install the [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) environment, the plotting programm [gnuplot](http://gnuplot.info/) as well as the IDE [geany](https://www.geany.org/) to your computer. In addition, it is going to compile and store the toolbox on your computer.


#### The installation process in detail

The file `install_fortran.bat` is a so-called batch-file, which we use to execute certain installation routines on your computer. You can run this file by double-clicking on it. During the installation process, we will make some changes to your system’s registry and (possibly) the file association register. For this, the file `install_fortran.bat` needs administrator rights, which you will have to grant by clicking `Yes` in the dialogue that shows up immediately when you execute the batch-file. Once administrator rights are granted, a console appears and you will be asked to confirm that you want to install all software components as well as our toolbox on your computer. You can do this by typing the letter `y` and pressing ENTER. We strongly recommend that you use the standard installation directory, which is `C:\cygwin`. If for some reasons you are not able to install Fortran to this directory, you can follow the advice in bullet point 3.

Now the installation process has started. In detail you will see the following things happening:

1. The first part of the installer copies [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) components to your system using [cywin](https://www.cygwin.com).
2. Next, we install [gnuplot](http://gnuplot.info/).
3. The third installer component will bring the IDE [geany](https://www.geany.org/) to your system.

When the three software components are installed, the installation batch-file will do a couple more things:

1. It will patch the software geany such that it already contains the right configuration for using Fortran with our toolbox.
2. It will compile and install the toolbox to a directory to which the compiler is pointed.
3. It will ask you whether all Fortran code files should be associated with our programming IDE geany. If you want this to happen, then just confirm by typing the letter y and pressing Enter.

**This is it!** If you saw no error messages appearing, then Fortran as well as our toolbox were successfully installed on your system.


#### Changing the installation directory

We strongly recommend that you use the standard installation directory, which is `C:\cygwin`. If for some reasons you are not able to install Fortran to this directory, you can change the installation directory according to the following steps:

1. You have to edit the file `install_fortran.bat`. To do so, right-click on the file and choose Edit from the drop-down menu that appears.
2. A text editor should appear that shows the content of the batch file. Go to line 12 of the batch file, which reads
    ```
    SET "location=C:\cygwin"
    ```
3. Now you can change `C:\cygwin` to the directory you want Fortran to be installed to. For example, if you want Fortran to be located in `D:\Fortran\myfortran`, you should change line 14 to
    ```
    SET "location=D:\Fortran\myfortran"
    ```
4. Save the file `install_fortran.bat` and exit the text editor.
5. Now run `install_fortran.bat` again. The installation directory should have changed.



## Uninstallation

In the unlikely event that you want to uninstall all the software components as well as the toolbox from your system, you can use the file `uninstall_fortran.bat`. It follows the same logic as the installation batch file. You just have to follow the instructions. Don’t forget to change the location of your Fortran installation, if you have changed the installation directory according to bullet point 3.


## Toolbox Update

Sometimes you might feel that your toolbox is out of date. In this case, just pull the most recent version of this github repository and execute the file `update_toolbox.bat`. This will compile the latest toolbox version and store is such that you can use it in all your programs.


## Run your first Fortran program

If you've successfully installed Fortran on you computer, [you can learn here how to run your first program](https://www.ce-fortran.com/run-program-win/).
