TRACS: Transient Current Simulator
===================================

Simulator of transient currents and charge collection in semiconductor detectors based on Shockley–Ramo theorem

![Screenshot of ](/docs/images/tracs_potentials_screenshot.png?raw=true)


# Requisites

  - C++11 compatible compiler
  - CMake
  - [Fenics Libraries](http://fenicsproject.org/download/)
  - Qt4
  - ROOT framework (for exporting histograms and generation of arbitrary charge distributions)

# Installation

    git clone https://github.com/pablodecm/TRACS
    cd TRACS
    mkdir build
    cmake ..
    make

