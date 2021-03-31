# Docker repository: [ce-fortran](https://hub.docker.com/r/fabiankindermann/ce-fortran)

This docker repository contains a basic version of the ce-fortran compiler toolkit. It is based on Ubuntu 20.04. In particular, it features the following packages:

1. The [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) for compiling Fortran source code.
2. [GNU GDB](https://gnu.org/software/gdb/) for debugging Fortran programs.
3. [GNU Make](https://www.gnu.org/software/make/) for creating more complicated compilation processes.
4. [gnuplot](http://gnuplot.info/) for creating graphs from Fortran programs. Note that the gnuplot installation in this repository does not have a graphical user interface, but only creates graphs in pdf, eps or png format.
5. The ce-fortran toolbox.


## Installation

In order to install the `ce-fortran` Docker repository, you need to have [Docker](https://www.docker.com) running on your computer. Please refer to the [Docker website](https://www.docker.com) for details on the installation process.

## Running the ce-fortran-vnc repository

In order to run this repository, open a terminal (or `cmd.exe` on Windows) and type
```Docker
docker run --rm -it fabiankindermann/ce-fortran
```
The base directory in this image is `/home/user`. It contains the source codes to both our [book](https://global.oup.com/academic/product/introduction-to-computational-economics-using-fortran-9780198804406?prevSortField=1&sortField=8&start=0&resultsPerPage=20&prevNumResPerPage=20&lang=en&cc=no) as well as its solution manual. You can start playing around with compiling and running source code in this repository. For example, you can navigate to the folder
```bash
cd /home/user/code-book/prog01/prog01_01
```
If you now want to run the program `prog01_01.f90`, simply type
```bash
gfortran $FCOMPILE -o prog prog01_01.f90 && ./prog
```
You can run all programs like this. `$FCOMPILE` is a preset environment variable that contains all necessary compiler flags and ensures that the toolbox is available.

## Accessing content from your own computer

You might want to start compiling source codes that are stored on your computer instead of the container. To this end, you will have to mount a local directory into the Docker container. This is a fairly simple task. 

**If you are on Linux or on macOS:**
Use the terminal to navigate to the folder in which your source code is stored. Then run
```Docker
docker run --rm -it -v $(pwd):/code -w /code fabiankindermann/ce-fortran
```


**If you are on Windows:**
Use the command prompt (`cmd.exe`) to navigate to the folder in which you source code is stored. Then run
```Docker
docker run --rm -it -v %cd%:/code -w /code fabiankindermann/ce-fortran
```

You will now enter the container in the directory `code`, which contains the file and subfolders of your current directory. You can execute any Fortran source code (or makefile) as described above.

## Alternative ways of using this container

There are ample other ways of using this container. Please refer to the [Docker website](https://www.docker.com) and the many tutorials that are out there to explore all features of Docker containers.`