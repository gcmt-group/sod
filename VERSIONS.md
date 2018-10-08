## Version 0.44 (October 2018)

- Added a simple pair-based extrapolation (spbe) calculator to predict energies for n>2 substitutions from the energies for n=0, 1, and 2 substitutions. 
- Makefile added; just type ```make all``` from ROOTSOD to compile all the fortran code.
- All scripts converted to bash; they are all in the ROOTSOD/bin folder, and have extension .sh
- Resolved issue in scripts related to number of columns from ```ls -al``` command.
- ```sod_stat``` now calculates thermodynamic properties at T=1K, 300K, 1000K and in the limit of a very high temperature, if the TEMPERATURES file does not exist.
 
