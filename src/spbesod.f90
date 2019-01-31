!*******************************************************************************
!    Copyright (c) 2018 Ricardo Grau-Crespo, Said Hamad
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!******************************************************************************

    PROGRAM spbesod 
       IMPLICIT NONE


       INTEGER :: i,j,k,p,q,m,m1,m2,m1ref,m2ref,mmin,mmax,Mm,Mm1,Mm2,op,nop,npos,pos,nsubs,aux,irescale
       REAL (KIND=8) :: E0,EN,EN_1,mu1,mu2,E1ref,E2ref,emin,emax
       CHARACTER*20 :: auxstring
       LOGICAL:: INSPBE_exists
       REAL (KIND=8),DIMENSION(:),   ALLOCATABLE:: dE1
       REAL (KIND=8),DIMENSION(:,:), ALLOCATABLE:: dE2
       REAL (KIND=8),DIMENSION(:),   ALLOCATABLE:: energies1
       REAL (KIND=8),DIMENSION(:),   ALLOCATABLE:: energies2
       REAL (KIND=8),DIMENSION(:),   ALLOCATABLE:: a1, a2, energies
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: eqmatrix
       INTEGER,DIMENSION(:), ALLOCATABLE:: omega, omega1, omega2
       INTEGER,DIMENSION(:), ALLOCATABLE:: conf1
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: conf2
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: conf


!!!!!! Input files

       OPEN (UNIT=10,FILE="ENERGIES0")
       OPEN (UNIT=11,FILE="ENERGIES1")
       OPEN (UNIT=12,FILE="ENERGIES2")
       OPEN (UNIT=13,FILE="EQMATRIX")
       OPEN (UNIT=14,FILE="OUTSOD")
       OPEN (UNIT=15,FILE="OUTSOD1")
       OPEN (UNIT=16,FILE="OUTSOD2")

       INQUIRE(FILE="INSPBE", EXIST=INSPBE_exists)
       if (INSPBE_exists) then
          OPEN (UNIT=17, FILE="INSPBE")
       else
          OPEN (UNIT=18, FILE="INSPBE.tmp")
       endif

!!!!!! Output files

       OPEN (UNIT=50,FILE="ENERGIES")
       OPEN (UNIT=99,FILE="OUTSPBE")


! DEFINITION OF VARIABLES:
!
! nsubs               numbers of substitutions
! op                  Index for the operators in the supercell
! nop                 Total number of operators in the supercell
! npos		      Number of atoms of the target species 
! pos		      Index for the positions
! conf		      List of all independent configurations for the structure with nsubs substitutions 
! conf1		      List of all independent configurations for the structure with 1 substitutions 
! conf2		      List of all independent configurations for the structure with 2 substitutions 
!                     All configurations are given as a list of the substituted positions
! m		      Index for the configurations (conf)     
! Mn                  Total number of configurations in conf (m=1,Mm)
! m1		      Index for the configurations (conf1) for the structure with 1 defect    
! Mm1                 Total number of configurations in conf1 (m1=1,Mm1) for the structure with 1 defect    
! m2		      Index for the configurations (conf2) for the structure with 2 defects    
! Mm2                 Total number of configurations in conf2 (m2=1,Mm2) for the structure with 2 defects    
! mu1		      Scaling parameter for single defect energies 
! mu2		      Scaling parameter for pair-wise interactions
! EN                  Energy of the x=1 endmember 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the EQMATRIX file 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       read(13,*) nop, npos
       ALLOCATE(eqmatrix(1:nop,1:npos))

       Do op=1,nop
          read(13,*) eqmatrix(op,1:npos)
       enddo



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the OUTSOD files:
!              - OUTSOD0 is not read, since we know the energy is E0
!              - OUTSOD1 for the structure with 1   substitution
!              - OUTSOD2 for the structure with 2   substitutions
!              - OUTSOD  for the structure with n>2 substitutions
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!       OUTSOD
!
        read(14,*)  nsubs
        read(14,*)  Mm

        ALLOCATE(conf(1:Mm,1:nsubs))
        ALLOCATE(omega(1:Mm))

        do m=1,Mm
           read(14,*)  aux,omega(m),conf(m,1:nsubs)
        enddo


!       OUTSOD1
!
        read(15,*)  
        read(15,*) Mm1 

        ALLOCATE(conf1(1:Mm1))
        ALLOCATE(omega1(1:Mm1))

        do m=1,Mm1
           read(15,*)  aux,omega1(m),conf1(m)
        enddo

!       OUTSOD2
!
        read(16,*)  
        read(16,*)  Mm2

        ALLOCATE(conf2(1:Mm2,1:2))
        ALLOCATE(omega2(1:Mm2))

        do m=1,Mm2
           read(16,*)  aux,omega2(m),conf2(m,1:2)
        enddo



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the files with the energies of 0, 1 and 2 substitutions (ENERGIES0, ENERGIES1 and  ENERGIES2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       ALLOCATE(energies1(1:Mm1))
       ALLOCATE(energies2(1:Mm2))

       read (10,*)   E0

       do m=1,Mm1
          read (11,*) energies1(m)
       enddo

       do m=1,Mm2
          read (12,*) energies2(m)
       enddo

       write(99,*) "-----------------------------------------------------------------"
       write(99,*) "Parameters calculated from data for 0, 1 and 2 substitutions:"
       write(99,*) "-----------------------------------------------------------------"
       write(99,*)
       write(99,*) "E0 (eV)"
       write(99,*) E0 



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating dE1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       ALLOCATE(dE1(1:npos))

       Do m1=1,Mm1
          Do op=1,nop
             dE1(eqmatrix(op,conf1(m1)))=energies1(m1)-E0
          Enddo
       Enddo


       write(99,*)
       write(99,*) "dE1 (eV)"
       write(99,*) dE1(:) 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating dE2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       ALLOCATE(dE2(1:npos,1:npos))

       Do m2=1,Mm2
          Do op=1,nop
             dE2(eqmatrix(op,conf2(m2,1)),eqmatrix(op,conf2(m2,2)))= &
             energies2(m2)-E0-dE1(eqmatrix(op,conf2(m2,1)))-dE1(eqmatrix(op,conf2(m2,2)))
          Enddo
       Enddo

       write(99,*)
       write(99,*) "dE2 (eV)"
       do p=1,npos
          write(99,*) dE2(p,1:npos) 
       enddo
       write(99,*)




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating configuration-dependent summations 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       ALLOCATE(a1(1:Mm))
       ALLOCATE(a2(1:Mm))

       Do m=1,Mm
          ! Initializing
            a1(m)=0.0
            a2(m)=0.0

          Do p=1,nsubs
             a1(m)=a1(m)+dE1(conf(m,p))
          Enddo

          Do p=1,nsubs-1
             Do q=p+1,nsubs
                a2(m)=a2(m)+dE2(conf(m,p),conf(m,q))
             Enddo
          Enddo

       Enddo
       
       write(99,*)


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSPBE file, if it exists, and calculating rescaling parameters mu1 and mu2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       if (INSPBE_exists) then

          read (17,*)
          read (17,*) irescale
          read (17,*)
          select case (irescale)
               case (0)
                    read (17,*)
                    write(*,*) "No rescaling (mu1=mu2=1.0)."
                    write(99,*) "No rescaling (mu1=mu2=1.0)."
                    mu1=1.0
                    mu2=1.0
               case (1)
                    read (17,*) mu1, mu2
                    write(*,*) "Rescaling parameters mu1 and mu2 read from INSPBE."
                    write(99,*) "Rescaling parameters mu1 and mu2 read from INSPBE."
                    write(99,*) mu1, mu2
               case (2)
                    read (17,*) m1ref, E1ref
                    read (17,*) m2ref, E2ref
                    mu1=(a2(m2ref)*(E1ref-E0)-a2(m1ref)*(E2ref-E0))/(a1(m1ref)*a2(m2ref)-a2(m1ref)*a1(m2ref))
                    mu2=(a1(m1ref)*(E2ref-E0)-a1(m2ref)*(E1ref-E0))/(a1(m1ref)*a2(m2ref)-a2(m1ref)*a1(m2ref))
                    write(*,*) "Rescaling parameters calculated using reference energies for two configurations."
                    write(99,*) "mu1 = ", mu1, "mu2 = ", mu2
               case default
                    mu1 = 0.0
                    mu2 = 0.0
                    write(*,*) "Invalid case for irescale in INSPBE - no rescaling will be applied (mu1=mu2=1.0)."
          end select

       else

          mu1 = 1.0
          mu2 = 1.0
          write(*,*) "No INSPBE file, therefore no rescaling applied (mu1=mu2=1.0)."
          write(99,*) "No INSPBE file, therefore no rescaling applied (mu1=mu2=1.0)."

       endif



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculating ENERGIES
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       ALLOCATE(energies(1:Mm))


          write(99,*)
          write(99,*)
          write(99,*) "Zeroth-, first-, and second-order contribution by configuration for ", nsubs, "substitutions:"
          write(99,*)
          write(99,*) "   m     E0(eV)             dE1(eV)         dE2(eV)      E[total] (eV)"
          write(99,*) " --------------------------------------------------------------------------"

       Do m=1,Mm
          energies(m)= E0 + mu1*a1(m) + mu2*a2(m)
          write(99,'(1I5,4F16.6)') m, E0, mu1*a1(m), mu2*a2(m),  energies(m)
          write(50,*) energies(m)
       Enddo

       write(99,*)
       write(99,*)

       mmin=MINLOC(energies(1:Mm),DIM=1)
       emin=MINVAL(energies(1:Mm))
       mmax=MAXLOC(energies(1:Mm),DIM=1)
       emax=MAXVAL(energies(1:Mm))

       write(99,*) "Minimum-energy configuration: ", mmin, " with energy: ", emin, " eV."
       write(99,*) "Maximum-energy configuration: ", mmax, " with energy: ", emax, " eV."
       write(99,*)
       write(99,*)

     if (.NOT.(INSPBE_exists)) then
      write(99,*) "To improve agreement with reference set, calculate reference energies for these two configurations," 
      write(99,*) "then input reference energies in INSPBE file. Use INSPBE.tmp as template." 
      write(99,*)
      write(18,*) "# irescale case: 0 = no rescaling; 1 = enter mu1 and mu2 manually; 2= enter two reference energies &
                  &(e.g. from DFT)"
      write(18,*) "2" 
      write(18,*) "# If irescale=1, enter one line with mu1, mu2; if irescale=2, enter two lines (m1, E1), (m2, E2)"
      write(18,*) mmin, emin 
      write(18,*) mmax, emax 
     endif


!!!!!!Deallocating arrays
       DEALLOCATE(eqmatrix)
       DEALLOCATE(conf)
       DEALLOCATE(omega)
       DEALLOCATE(conf1)
       DEALLOCATE(omega1)
       DEALLOCATE(conf2)
       DEALLOCATE(omega2)
       DEALLOCATE(energies1)
       DEALLOCATE(energies2)
       DEALLOCATE(dE1)
       DEALLOCATE(dE2)
       DEALLOCATE(energies)
       DEALLOCATE(a1)
       DEALLOCATE(a2)


!!!!!!Reporting the end
       write (*,*) "Done!!!"
       write (*,*) 
       write (*,*) 

    END PROGRAM spbesod

