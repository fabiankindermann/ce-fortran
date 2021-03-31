# Docker repository: [ce-fortran-vnc](https://hub.docker.com/r/fabiankindermann/ce-fortran-vnc)

This docker repository contains the ce-fortran compiler toolkit in a Ubuntu 20.04 container with VNC access. The repository is built on [dorowu/ubuntu-desktop-lxde-vnc](https://hub.docker.com/r/dorowu/ubuntu-desktop-lxde-vnc/). It features the following packages:

1. The [GNU Fortran Compiler](https://gcc.gnu.org/fortran/) for compiling Fortran source code.
2. [GNU GDB](https://gnu.org/software/gdb/) for debugging Fortran programs.
3. [GNU Make](https://www.gnu.org/software/make/) for creating more complicated compilation processes.
4. [gnuplot](http://gnuplot.info/) for creating graphs from Fortran programs with an X11 terminal support.
5. The IDE [geany](https://www.geany.org/), which is patched to work with the gfortran installation as well as the ce-fortran toolbox.
6. The ce-fortran toolbox.

## Installation

In order to install the `ce-fortran-vnc` Docker repository, you need to have [Docker](https://www.docker.com) running on your computer. Please refer to the [Docker website](https://www.docker.com) for details on the installation process.

## Running the ce-fortran repository

In order to run this repository, open a terminal (or `cmd.exe` on Windows) and type
```Docker
docker run --rm -it -p 6080:80 fabiankindermann/ce-fortran-vnc
```
If the docker cotainers has started successfully, use a browser with VNC support and open [127.0.0.1:6080](). You should now see the graphical user interface of the Ubuntu system running in your container. The base directory in this container is `/root`. It contains the source codes to both our [book](https://global.oup.com/academic/product/introduction-to-computational-economics-using-fortran-9780198804406?prevSortField=1&sortField=8&start=0&resultsPerPage=20&prevNumResPerPage=20&lang=en&cc=no) as well as its solution manual. If you've successfully logged into the container, [you can learn here how to run your first program](https://www.ce-fortran.com/run-program-lin/).

## Accessing content from your own computer

You might want to start compiling source codes that are stored on your computer instead of the container. To this end, you will have to mount a local directory into the Docker container. This is a fairly simple task. 

**If you are on Linux or on macOS:**
Use the terminal to navigate to the folder in which your source code is stored. Then run
```Docker
docker run --rm -it -p 6080:80 -v $(pwd):/root/code fabiankindermann/ce-fortran-vnc
```

**If you are on Windows:**
Use the command prompt (`cmd.exe`) to navigate to the folder in which you source code is stored. Then run
```Docker
docker run --rm -it -p 6080:80 -v %cd%:/root/code fabiankindermann/ce-fortran-vnc
```

There will now be a directory `code` in the containers home folder, which contains the file and subfolders of your current directory.

## Alternative ways of using this container

There are ample other ways of using this container. Please refer to the [Docker website](https://www.docker.com) and the many tutorials that are out there to explore all features of Docker containers.
