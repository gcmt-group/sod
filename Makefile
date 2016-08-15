## This makefile must be executed with gmake (gnu make).
all:
	gfortran -o bin/statsod src/statsod/statsod.f90
	gfortran -o bin/combsod   src/combsod/bubble.f90  src/combsod/ksubset.f90  src/combsod/member.f90  src/combsod/cell.f90 src/combsod/ccf.f90 src/combsod/combsod.f90

clean:
	rm bin/* 

