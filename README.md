*******************************************************************************
                SOD Version 0.40 - Notes for users
*******************************************************************************

SOD (standing for Site-Occupation Disorder) is a package of tools for the computer modelling of periodic systems with site disorder, using the supercell ensemble method. 

The package is distributed under the GPL licence. Please see details in the COPYING file. 

You can find below the essential info needed to use SOD . Please note that I can give only limited support to users.

1. Functionalities:

- identification of all inequivalent configurations of site substitutions in an arbitrary supercell of an  initial target structure with any space group symmetry.

- calculation of the degeneracies of configurations

- generation of input files for codes like Gulp and Vasp

- statistical mechanics processing of output

2. Installing SOD:

- Download the file sod(version).tar (e.g. sod0.26.tar) and copy to a directory, say ROOTSOD

- tar xvf sod(version).tar

- add ROOTSOD/sod(version)/exe to your executables path (e.g. in bash, add the line 'export PATH=$PATH:ROOTSOD/sod(version)/exe' to your .bashrc file)

3. Content of the folders:

- sod(version)/src contains the source files.
- sod(version)/sgo is a library of space group operators (e.g. 131.sgo contains the operators of the space group 131).
- sod(version)/exe contains the executables. Linux executables are provided here.
- sod(version)/examples contains three examples, based on the cubic perovskite, rutile and rocksalt structures.

4. Running SOD

- We recommend to create a new folder (say FOLDER_NAME) for each application. This will be referred to as the working directory.

- In FOLDER_NAME, you must create a file named INSOD which contains all the information for running the combinatorics part of the program. Use the INSOD file given in one of the examples. The file is self-explanatory. In the current version the format of this file is rigid, so keep the same number of blank lines.

- In FOLDER_NAME, you must also include a file named SGO with the matrix-vector representations of the symmetry operators. First check if your space group is included in the ROOTSOD/sod(version)/sgo library; if this is the case, just copy the file into your working directory, under the name SGO:

cp ROOTSOD/sod(version)/sgo ./SGO

otherwise you have to create the file using the Tables of Crystallography, or from the website of the Bilbao Crystallographic Server <www.cryst.ehy.es>. The first three numbers in each line are one row of the operator matrix and the fourth number is the component of the operator's translation vector.

- If you want to generate Gulp input files for all the independent configurations found by SOD, in addition to setting FILER=1 in the INPUT file, you must provide two files in the working directory:

top.gulp contains the heading of the gulp input file (until the keyword cell).

bottom.gulp contains the tail of the gulp input file (everything after the list of coordinates, including species, potentials, etc).

- To run the combinatorics program, just type:

sod_comb

5. Output of the sod_comb programme.

- When the programme finishes, it writes to the standard output the total number of configurations and the number of independent configurations according to the crystal symmetry, plus some other useful information.

- It also writes the data file OUTSOD, which contains information on the independent configurations (one line for each configuration). The first number is the index of the configuration, the second is its degeneracy, and the next numbers are the substitution sites.

- The directory CALCS is generated, which contains the input files for Gulp or VASP  and a script that sends the job.

6. Configurational averages and thermodynamics:

In order to do configurational averages and thermodynamics, you need to execute the program

sod_stat

which requires 4 input files:

- OUTSOD, which was the output from sod_comb

- TEMPERATURES, a list of temperatures for the Boltzmann statistics, in one column, e.g.:

600
1000
1500
2000

- ENERGIES, which contains (in one column) the energies of all the configurations, in the same order that they were generated by SOD (like in the OUTSOD file). There are some scripts in ROOTSOD/sod(version)/exe/  to help you do this:

a) If you are using GULP, the script sod_gulp_ener will extract all the energies, assuming all output files,  with extension .gout, are in the same folder. If you have calculated vibrational free energies for each configurations, sod_gulp_free will extract these. Depending on your Unix system, you might need to change $9 (for example, to $8) in sod_gulp_ener (this refers to the number of the column containing the file names when you execute the command "ls -al" in your unix system).

b) If you are using VASP, the script sod_vasp_ener will extract all the energies, assuming you have folders 01, 02, 03, etc. each containing one calculation.  If you have more than 99 configurations and therefore have folders 001, 002, ..., 131, 132, etc. then you should edit the sod_vasp_ener script, changing ?? to ???. Depending on your Unix system, you might need to change $9 (for example, to $8) in sod_gulp_ener (this refers to the number of the column containing the file names when you execute the command "ls -al" in your unix system).

- DATA, which contains n colums of data to average. The first line contains just the number n of columns to read. For example:

2
34.5   4.34
37.7   4.35
35.6   4.38
38.8   4.41

The data can be cell lenghts, or volumes (please see SOD papers for strategies on how to obtain average cell parameters) or any other observable obtained from the calculations. Scripts like sod_vasp_cell can help you do this, please edit carefully before using them.

sod_stat will generate two files: probabilities.dat and statistics.dat, whose content is self-explanatory.

Very important note: While configurational averages (e.g. of cell parameters and enthalpies) tend to converge very quickly with supercell size, entropies and free energies, which are not defined by averaging, converge very slowly with supercell size, and are generally in large error when using the SOD method. I therefore do not recommend using SOD for the calculation of entropies and free energies, unless appropriate correcting procedures have been used.



7. If you use SOD in your research work, please include a citation to this article:

"Symmetry-adapted configurational modelling of fractional site
occupancy in solids." R. Grau-Crespo, S. Hamad, C.R.A. Catlow and N. H. de Leeuw. Journal of
Physics - Condensed Matter 19, 256201 (2007)

and if possible send me a pdf copy of your paper.


Happy SODing!!!

Ricardo Grau-Crespo

