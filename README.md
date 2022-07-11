<img src="docs/HFCalc.png"></img>

<b>HFCalc</b> is a first-principles software package that can currently calculate total energy of molecules using restricted Hartree-Fock theory. It is written in C++. It employs linear combination of atomic orbitals, where the atomic orbitals are represented using contracted Gaussians (STO-NG basis set). It is free, open-source and in constant development. Header files from the Eigen and Boost libraries are used.

## Installation

The package can be installed using make.
```
cd src/
make CXX=... CXXFLAGS=...
make install
```
The executable (called hfcalc) will be placed in the bin/ folder. Note that the boost header files give an error when using Apple clang compiler if you don't specify -std=c++11 in the CXXFLAGS. You are welcome to try other compilers, please open an issue if it doesn't work.

## Running a calculation

To run a calculation, we need two files, an input file called config.in, and a position file the name of which has to be specified in the input file. Let's go through a calculation for the hydrogen molecule. The config.in file looks like this
```
POS posdata.dat
BASIS STO3G
  # Current options are STO2G,STO3G,STO4G,STO5G,STO6G  
BFOLDER ../../basis/
MODE 0
  #  0 - Restricted Hartree-Fock
  #  1 - Unrestricted Hartree-Fock (in progress)
RELAX 0
  #  0 - no relaxation
  #  1 - relax structure
MIXTYPE 0
  #  0 - Simple Mixing
  #  1 - Anderson Mixing (in progress)
  #  2 - Broyden Mixing (in progress)
WRITE 0
  #  0 - write out energies only (default)
  #  1 - write out wavevector and contracted gaussian coefficients
MIXPARAM 1
NITER 100
ETOL  0.000001
OUT outfile
```
Note that in this example, the current working directory is assumed to be sample/H2, hence BFOLDER (path to basis folder) is given as  ../../basis/. If you are running the calculation elsewhere, you must modify this input parameter accordingly. Also the '/' at the end of the path is important. <br>
The position file (posdata.dat) looks like that (positions are in atomic units). The second column is atomic number and third column is number of electrons. They are different in order to allow ions.
```
H 1 1 0.0 0.0 0.0  
H 1 1 0.0 0.0 1.4 
```
In the same directory, run the executable
```
path_to_repository/bin/hfcalc
```
The standard output will print the final converged energy
```
Execution started at Sun Jul  3 18:41:07 2022
###################################################################################

***       ***    ***********   **********          *          ***        **********
***       ***    ***********   **********         ***         ***        **********
***       ***    ***           **                *****        ***        **
*************    *******       **               *** ***       ***        **
*************    *******       **              ***   ***      ***        **
***       ***    ***           **             ***********     ***        **
***       ***    ***           **********    ***       ***    *********  **********
***       ***    ***           **********   ***         ***   *********  **********

###################################################################################
Version 1.0 (July 2022)
###################################################################################
Total Converged Energy : -1.11671 Hartrees
Data written to : outfile
```
The outfile contains more detailed information on the calculation. More examples are available in the sample folder. You are encouraged to run those examples and other systems of interest. Any unphysical or problematic results should be reported by opening an issue.

## Acknowledgement

This code was written by following the following
 * Thijssen, J. (2007). Computational Physics (2nd ed.). Cambridge: Cambridge University Press. doi:10.1017/CBO9781139171397
 * Szabo, A., & Ostlund, N. S. (1996). Modern quantum chemistry: Introduction to advanced electronic structure theory. Mineola, N.Y: Dover Publications.
 * https://joshuagoings.com/assets/integrals.pdf (for molecular integrals)
