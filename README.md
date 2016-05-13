TRACS: TRAnsient Current Simulator
===================================

Simulator of transient currents and charge collection in semiconductor detectors based on Shockley–Ramo theorem

![Screenshot of ](/docs/images/TRACS_CNVT.png?raw=true)

# Requisites

  - [C++11 compatible compiler](http://en.cppreference.com/w/cpp/compiler_support) (**Important!**) ||
    Tested only with g++-4.8 and Clang 6.x

  - [CMake](http://www.cmake.org/download/)
  
  - [Fenics Libraries](http://fenicsproject.org/download/) ||
    Compatibility appears to be broken for 1.6.x and up
  
  - [Qt4](http://download.qt.io/archive/qt/) ||
    Qt5 is not compatible 

  - [ROOT](https://root.cern.ch/downloading-root) ||
    Used for exporting histograms and generation of arbitrary charge distributions
    Recommended version 6.x

# Installation

Once you have installed all the software listed under "Requisites" you can now begin the installation process:

1) Get the source

    git clone https://github.com/IFCA-HEP/TRACS
    
2)[OPTIONAL-Recommended] Recompile fenics files

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
    
5) Compile in a machine with [N] cores

    make -j[N]
    
6) Move necessary files to the directory from which you will execute TRACS. eg.: TRACS/build/bin

    mv ../src/files2move2bin/* ./bin/
   
7) Execute TRACS

    cd bin/
    ./TRACS # for Command line version
    ./TRACS-GUI # for Grafical User Interface version

# Brief Introduction on How TRACS Works

  A summary of how TRACS functions can be found here. As a brief introduction it is meant for young padowams that have never seen TRACS before. More experienced users or anyone with some experience in C++11 or detector simulation might find it more useful to read the code and comments in folder TRACS/src/                                
                                                                       
  If you decide to use/develop TRACS you will be part of a great scientific project but you should proceed at your own risk as the only thing TRACS usage/development guarantees is headeaches and despairation.                                                         

# User's Guide for TRACS-GUI

## Typical use case scenario
						   
A typical user of TRACS-GUI would very likely wish to do some quick simulations in order to check the validity of certain parameters (e.g. Neff, trapping time...). For this checks one might want to take a look at the electric and weighting fields, do some sample transient current simulations and maybe even do some actual laser illumination  simulation. After the simulations are done one might want to change the aforementioned parameters and redo the simulations. Here we present a step-by-step guide for the user to achieve said checks succesfully without any previous knowledge required.
         
  Though it is not compulsory to follow the step in the apropiate order, it is recomended at least the first time to follow this process as explained for it helps to progresive checks where the next checks depends primarily on the previous item checked been correct and simulations/plots grow in complexity as we one from one step to another and so does the computing time associated with said simulations.



### 1) Setting up the detector
  ![Screenshot of ](/docs/images/TRACS_CNVT.png?raw=true)

  The detector is the main component of TRACS simulation around which the rest of the items pivot. 

 To properly set up the detector one should go to the potentials tab. In this tab electric and weighting potentials will be plotted. It also holds the input boxes to set all the properties of the detector. This features include size of the detector and it's strips, temperature, doping types and depletion voltage. In this tab one can also change the FEM solver properties such as bias voltage applied to the detector, mesh size (number of cells) and even the number of neighbouring strips [nns] in the microstrip one wishes to simulate (set to 0 for diode simulation). 

 Once the correct values have been inserted, one should click on the "Solve FEM Problem" button to initialize the detector and calculate its fields and potentials. From this tab one can see 2D and 3D prots of the electric and wheighting potentials.

### 2) Cheking the fields
 ![Screenshot of ](/docs/images/TRACS_Fields_DJ.png?raw=true)


  After the detector has been initialized and FEM problem has been solved one can check the shape of the fields in the fields tab. Here the user is greeted with a 2D pot of the fields on the left hand side buttons in the middle and empty plots for slicing the fields in a 1D plot on the right hand side.


 The user can decide to plot X, Y components or its modulus as well as to define the point at which the slice is taken. 2D and 3D plots ar also available throught VTK.

### 3) Fast transient checks
![Screenshot of ](/docs/images/TRACS_Curr_DP.png?raw=true)

  Next logical step is to check that the currents araising from the previously mentioned fields are as one expected. For that we move on to the currents tab.

  Here we have a 2D plot of the detector with a coloured superposition of the electric field inside of it. The fastest check that can be done is performed by double clicking anyonwhere inside the detector. This will create an electron-hole pair and will drift it throught the detector. The current given by this pair will be ploted below the detector showing separated and joint contributions. One can also change the total charge simulated as well as numerically imput the position of the pair. This is an inmediatecheck (simulation) but very limited in term of real-world comparison. 

 For a more real-world scenario that is still almost inmediate one can use the line drifting tool below. The user can define the starting and ending points of a line in which the carriers will be placed and drifted.
 
 This allows for more complex simulations such as edge-TCT configuration (a horizontal line) or red-TCT (short segment near the ends of the detector) whilst still being almost inmediate to simulate. The total amount of charge drifted as well as the number of carriers drifted can be modified by setting the separation between carriers in the line. A typical value around "Separation = 10" has been found to give almost inmediate results while maintaining great accuracy.

  An common option for simulating the electronics shaping is available for both single pair and line simulations. The input parameter is the capacitance of the Resistor-Capacitor circuit that is simulated using a simple "low pass filter" type of loop. The resistance is assumed to be 50 Ohm as is the case in almost all the real-world cases.

### 4) Simulating real-world scenarios
![Screenshot of ](/docs/images/TRACS_Carr_DP.png?raw=true)

  If all this tests are not enough, TRACS-GUI also provides a more sofisticated simulation including an arbitrary carrier distribution that can be made to simulate any laser illumination. Together with the source code TRACS provides various sample distributions including red-top-TCT and 3 heights of edge-TCT.

  For this kind of simulation the user should go to the Carriers tab. Here one can input the desired carriers file and load it to the program for the simulation. When loaded the carrier distribution will appear as a 2D plot on the left hand side of the window so that the user can check that the carrier distribution is loaded correctly before simulating.

  The user should then adjust the correct time to simulate and the time step to get the desired compromise between computing time and accuracy. Test show that a while a 50ps binning is recommended for full simulations, a time step of 100ps or even higher works good for GUI simulations while taking less computing time. The option for simulating the RC-shaping of the signal is also available with capacitance input that is independent of the one mentiones in step three. 

### 5) Changing simulation parameters

  This last simulation reflects the real-world case that would be simulated using the CLI interface. If the user is content with the results he/she should note the results down and modify the "Config.TRACS" file to run a more complete simulation using the CLI version.

  If the results are not the desired ones, the user should now modify the parameters and check if the new values yeild better results. For succesfully making changes to the original values there are some necessary steps that might not appear evident to the new user; we list them below: 

 - After modification of any of the detector parameters the user should click on the Solve FEM Problem again to load them. This is due to the detector only been updated  when the button is clicked and FEM problem solved again.

- Carriers file should be (re)loaded before every simulation in the carriers tab. This ensures that the the detector used for the simulation is updated to the one set by solving the FEM Problem. 

  This process can be repeted as many times as desired though we advise to close and rerun the program should any problems appear. A known issue is that in the Carriers tab it some times takes various attempts for the 
simulation to reflect the new configuration after modifying parameters; this problem is solved closing and rerunning TRACS-GUI.
