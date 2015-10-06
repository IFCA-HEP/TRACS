TRACS: Transient Current Simulator
===================================

Simulator of transient currents and charge collection in semiconductor detectors based on Shockley–Ramo theorem

![Screenshot of ](/docs/images/TRACS_CNVT.png?raw=true)

# Requisites

  - [C++11 compatible compiler](http://en.cppreference.com/w/cpp/compiler_support) (**Important!**) ||
    Tested only with g++-4.8

  - [CMake](http://www.cmake.org/download/)
  
  - [Fenics Libraries](http://fenicsproject.org/download/)
  
  - [Qt4](http://download.qt.io/archive/qt/) ||
    Qt5 is not compatible 

  - [ROOT](https://root.cern.ch/downloading-root) ||
    Used for exporting histograms and generation of arbitrary charge distributions
    Recommended version 6.x

# Installation

1) Get the source

    git clone https://github.com/IFCA-HEP/TRACS
    
2)[OPTIONAL] Recompile fenics files

    cd TRACS/src
    ffc -l dolfin Poisson.ufl
    ffc -l dolfin Gradient.ufl
    cd ../
    
3) Creat folder to store executables

    cd TRACS
    mkdir build
    cd build
    
4) Configure CMake and creat Makefile

    cmake ..
    
5) Compile in a machin with [N] cores

    make -j[N]
    
6) Move necessary files to the directory from which you will execute TRACS. eg.: TRACS/build/bin

    mv ../src/files2move2bin/* ./bin/
   
7) Execute TRACS

    cd bin/
    ./TRACS # for Command line version
    ./TRACS-GUI # for Grafical User Interface version

# Brief Introduction on How TRACS Works

![Screenshot of ](/docs/images/TRACS_CNVT.png?raw=true)
![Screenshot of ](/docs/images/TRACS_Fields_DJ.png?raw=true)
![Screenshot of ](/docs/images/TRACS_Carr_DP.png?raw=true)
![Screenshot of ](/docs/images/TRACS_Curr_DP.png?raw=true)


  A summary of how TRACS functions can be found here. As a brief introduction it is meant for young padowams that have never seen TRACS before. More experienced users or anyone with some experience in C++11 or detector simulation might find it more useful to read the code and comments in folder TRACS/src/                                
                                                                       
  If you decide to use/develop TRACS you will be part of a great scientific project but you should proceed at your own risk as the only thing TRACS usage/development guarantees is headeaches and despairation.                                                         
