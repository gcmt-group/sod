## This makefile must be executed with gmake (gnu make).

f90comp = gfortran

all:
	$(f90comp) -o bin/combsod src/factorials.f90 src/bubble.f90  src/ksubset.f90  src/member.f90  src/cell.f90 src/ccf.f90 src/combsod.f90
	$(f90comp) -o bin/genersod src/member.f90 src/cell.f90 src/genersod.f90
	$(f90comp) -o bin/spbesod src/spbesod.f90
	$(f90comp) -o bin/invertOUTSOD src/invertOUTSOD.f90
	$(f90comp) -o bin/statsod  src/statsod.f90
	$(f90comp) -o bin/gcstatsod  src/factorials.f90 src/momenta.f90 src/gcstatsod.f90
	$(f90comp) -o bin/peaks2spec  src/peaks2spec.f90

clean:
	rm bin/combsod  
	rm bin/invertOUTSOD  
	rm bin/statsod
	rm bin/gcstatsod
	rm bin/genersod
	rm bin/spbesod


