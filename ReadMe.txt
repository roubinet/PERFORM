//////////////////////
// Code description //
//////////////////////
// The code PERFORM is dedicated to simulate solute transport in heterogeneous fractured porous media.

// Simulation conditions are described in the folder "Input" where the file "File_names.txt" contains the names of the parameter file (located in the folder "Param_files") and the domain file (located in the folder "Domain_files"). 

// Examples of these files are present in each of these two folders corresponding to the simulations presenting in the article [Roubinet et al., 2012].

// The code requires the following libraries:
- CGAL-3.2.1
- boost_1_49_0

// For more information, please see the related documentation file and article [Roubinet et al., 2012] 


//////////////////////
// Code compilation //
//////////////////////
- create a new C++ project 
- import the source files located in the folder "src"
- add the following include folders 
(in Properties -> C/C++ General -> Paths and Symbols -> Includes):
your_path_to_libraries/CGAL-3.2.1/include
your_path_to_libraries/CGAL-3.2.1/include/CGAL/config/msvc7
boost_1_49_0
- add the library in libCGAL.a
(in Properties -> C/C++ General -> Paths and Symbols -> Libraries)

- .vcproj and .sln files for Visual Studio are located in the folder VStudio 

////////////////////
// Code execution //
////////////////////
// When running the code, user should enter the path containing the folders "Input" and "Output"

// Simulation results corresponding to solute arrival times will be located in the folder "Output" 


////////////////
// References //
////////////////
// D. Roubinet, J.-R. de Dreuzy, and D. M. Tartakovsky, Particle-tracking simulations of anomalous transport in hierarchically fractured rocks, Computer and Geosciences, 2012