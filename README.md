TRACS: Transient Current Simulator
===================================

Simulator of transient currents and charge collection in semiconductor detectors based on Shockley–Ramo theorem

![Screenshot of ](/docs/images/tracs_potentials_screenshot.png?raw=true)


# Requisites

  - C++11 compatible compiler (**Important!**)
  - CMake
  - [Fenics Libraries](http://fenicsproject.org/download/)
  - Qt4
  - ROOT framework (for exporting histograms and generation of arbitrary charge distributions)

# Installation

    git clone https://github.com/IFCA-HEP/TRACS
    cd TRACS
    mkdir build
    cd build
    cmake ..
    make
# How TRACS Works

===========================================================================
||                                                                       ||
||                BRIEF INTRODUCTION TO HOW TRACS WORKS                  ||
||                                                                       ||
||     A summary of how TRACS functions can be found here. As a brief    ||
|| introduction it is meant for young padowams that have never seen      ||
|| TRACS before. More experienced users or anyone with some experience   ||
|| in C++11 or detector simulation might find it more useful to read the ||
|| code and comments in folder TRACS/src/                                ||
||                                                                       ||
||    If you decide to use/develop TRACS you will be part of a great     ||
|| scientific project but you should proceed at your own risk as the     ||
|| only thing TRACS usage/development guarantees is headeaches and       ||
|| despairation.                                                         ||
||                                                                       ||
===========================================================================


=----------------------------- THE INPUTS --------------------------------=


=--------------------- SETTING UP THE SIMULATION--------------------------=

   
=--------------------------- SIMULATION --------------------------------=


=--------------------------- THE OUTPUT--------------------------------=
