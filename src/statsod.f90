!*******************************************************************************
!    Copyright (c) 2014 Ricardo Grau-Crespo, Said Hamad
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

       PROGRAM stats 

       IMPLICIT NONE

       INTEGER,PARAMETER :: NCONFMAX=10000, NCOLMAX=10, NPOINTSMAX=800, NTEMPMAX=50
       REAL (kind=8),PARAMETER :: kB=8.61734E-5, tolprob=1.0E-12, tolminspec=1.0E-6
       INTEGER :: m, auxm, Mm, g, ncol, npoints, col, tt, Ntt, nsubs, point
       INTEGER :: mEmin, mEmax
       REAL (kind=8)  :: Emin, Emax, paux, maxspec
       REAL (kind=8),DIMENSION(NTEMPMAX) :: Z,E,F,S,T
       REAL (kind=8) :: Zinf, Einf, Sinf
       REAL (kind=8),DIMENSION(NCONFMAX) :: ene, Erel 
       REAL (kind=8),DIMENSION(NCONFMAX,NTEMPMAX) :: p
       REAL (kind=8),DIMENSION(NCONFMAX) :: pinf
       INTEGER,DIMENSION(NCONFMAX) :: conf,omega
       REAL (kind=8),DIMENSION(NCOLMAX,NCONFMAX) :: data 
       REAL (kind=8),DIMENSION(NCOLMAX,NTEMPMAX) :: avedata 
       REAL (kind=8),DIMENSION(NPOINTSMAX,NCONFMAX) :: spec 
       REAL (kind=8),DIMENSION(NPOINTSMAX,NTEMPMAX) :: avespec 
       REAL (kind=8),DIMENSION(NCOLMAX) :: avedatainf 
       REAL (kind=8),DIMENSION(NPOINTSMAX) :: xspec, avespecinf 
       CHARACTER(len=30) fmtemplist
       LOGICAL :: TEMPERATURES_exists, DATA_exists, SPECTRA_exists
       INTEGER :: n, npos



!Input files

       INQUIRE(FILE="TEMPERATURES", EXIST=TEMPERATURES_exists)
       if (TEMPERATURES_exists) then
          OPEN (UNIT=10, FILE="TEMPERATURES")
       endif

       OPEN (UNIT=11,FILE="OUTSOD")
       OPEN (UNIT=12,FILE="ENERGIES")

       INQUIRE(FILE="DATA", EXIST=DATA_exists)
       if (DATA_exists) then
          write (*,*) "DATA file found"
          write (*,*) 
          OPEN (UNIT=13, FILE="DATA")
       else
          write(*,*) "DATA file not found. No averaging of DATA will be performed."
       endif

       INQUIRE(FILE="SPECTRA", EXIST=SPECTRA_exists)
       if (SPECTRA_exists) then
          write (*,*) "SPECTRA file found"
          write (*,*) 
          OPEN (UNIT=14, FILE="SPECTRA")
          OPEN (UNIT=15, FILE="XSPEC")
       else
          write(*,*) "SPECTRA file not found. No averaging of SPECTRA will be performed."
       endif

!Output files
       OPEN (UNIT=20,FILE="probabilities.dat")
       OPEN (UNIT=21,FILE="thermodynamics.dat")
       
       if (DATA_exists) then
           OPEN (UNIT=22,FILE="ave_data.dat")
       endif

       if (SPECTRA_exists) then
           OPEN (UNIT=23,FILE="ave_spectra.dat")
       endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Read the TEMPERATURES files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if (TEMPERATURES_exists) then
           tt=1
  1        continue
           read(10,*,end=10)  T(tt)
           tt=tt+1
           goto 1
  10       continue
           Ntt=tt-1
           CLOSE (10)
        else
           T(1)=1.0
           T(2)=300.0
           T(3)=1000.0
           Ntt=3
        endif



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read the OUTSOD file, which is the output from SOD, 
    !      giving configuration numbers and degeneracies
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    read(11,*) nsubs 
    read(11,*) Mm 
    do m=1, Mm
        read(11,*)  auxm, omega(m)
    enddo
    CLOSE (11)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read ENERGIES 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    do m=1,Mm
       read(12,*) ene(m)
    enddo
    CLOSE (12)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read the data file, DATA 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (DATA_exists) then

       read(13,*) ncol
       do m=1,Mm
          read(13,*) data(1:ncol,m) 
       enddo
       CLOSE (13)
    endif

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read the SPECTRA and XSPEC files, and calculate maximum spectrum intensity 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    if (SPECTRA_exists) then
       write(*,*) "Reading SPECTRA file..."
       write(*,*) 
       read(14,*) npoints
       do m=1,Mm
          read(14,*) spec(1:npoints,m)
       enddo
       CLOSE (14)


       maxspec=MAXVAL(spec)

       do point=1,npoints 
          read(15,*)  xspec(point)
       enddo
       CLOSE (15)

    endif


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       Start writing thermodynamics.dat file
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    write(21,*) "       T             E               F          S             "

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       Start writing ave_data.dat file
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (DATA_exists) then
       write(22,*) "       T    Average data"
    endif

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       Start writing ave_spectra.dat file
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (SPECTRA_exists) then
!       WRITE(fmtemplist,'(a4,i1,a15)') "(a6,",Ntt,adjustl("(f13.1,2x),a12)")
!       write(23,fmtemplist) "x   ",     T(1:Ntt),"    Infinity"
! XXX Arreglar esto
       write(23,*) "x   ",     T(1:Ntt)
    endif



            
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !          Get minimum and maximum energies and relative energies with
    !          respect to minimum in order to get better accuracy for the 
    !          exponential calculations
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    Emin=ene(1)
    Emax=ene(1)
    mEmin=1
    mEmax=1
    do m=2,Mm
       if (Emin > ene(m)) then
           Emin=ene(m)
           mEmin=m
       endif
       if (Emax < ene(m)) then
           Emax=ene(m)
           mEmax=m
       endif
    enddo                        

    write(20,*) "Configuration with minimum energy: ", mEmin
    write(20,*) "Configuration with maximum energy: ", mEmax
    write(20,*)

    do m=1,Mm
        Erel(m) = ene(m) - Emin
    enddo


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       This starts the loop over all temperature values
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    do tt=1,Ntt

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !          Calculate the partition function and probabilities
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Z(tt)=0.0
       do m=1,Mm
          Z(tt)=Z(tt)+omega(m)*exp(-Erel(m)/(kB*T(tt)))
       enddo
    
       do m=1,Mm
          p(m,tt)=omega(m)*exp(-Erel(m)/(kB*T(tt)))/Z(tt)
       enddo

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      Calculate the energy (E) and Helmholtz free energy (F)
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       E(tt)=0.0
       S(tt)=0.0
     
       do m=1,Mm
          E(tt)=E(tt)+ene(m)*p(m,tt)
       enddo

       F(tt)=Emin-kB*T(tt)*log(Z(tt))

       S(tt)=(E(tt)-F(tt))/T(tt)
    

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      Calculate the average data
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if (DATA_exists) then
        
           avedata(1:ncol,tt)=0.0
        
           do col=1,ncol
              do m=1,Mm
                 avedata(col,tt)=avedata(col,tt)+data(col,m)*p(m,tt)
              enddo
           enddo
        
        endif

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      Calculate the average spectra
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if (SPECTRA_exists) then

           avespec(1:npoints,tt)=0.0

           do point=1,npoints
              do m=1,Mm
                 avespec(point,tt)=avespec(point,tt)+spec(point,m)*p(m,tt)
              enddo
              if (avespec(point,tt)/maxspec.lt.tolminspec) then 
                 avespec(point,tt)=0.0
              endif
           enddo

        endif


       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !          Write PROBABILITIES for each temperature T 
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(20,100) "Temperature=",T(tt)
  100         format(a12,2x,f10.4)
       write(20,*) "        m  omega(m)   Erel(m)/eV       p(m)      p(m)/omega(m)"
       do m=1,Mm
          write(20,101) m,omega(m),Erel(m),p(m,tt),p(m,tt)/omega(m)
  101         format(i10,2x,i8,2x,e12.6,2(2x,f10.4))
       enddo

       write(20,*) 
       write(20,*)

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !          Write thermodynamic information for each temperature T 
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       write(21,200)  T(tt), E(tt), F(tt), S(tt)
  200          format(f10.1,2x,2(f14.4,2x),e12.6)



       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !       Writing ave_data for each temperature
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
       if (DATA_exists) then
          write(22,201)  T(tt), avedata(1:ncol,tt) 
  201          format(f10.1,2x,10(f10.4,2x))
       endif





    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       This ends the loop over all temperature values
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    enddo 

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       Calculating results in the limit of infinite temperature
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
    Zinf=0.0
    do m=1,Mm
       Zinf=Zinf+omega(m)
    enddo                                                                                
 
    do m=1,Mm
       pinf(m)=omega(m)/Zinf
    enddo
 
    do m=1,Mm
       Einf=Einf+ene(m)*pinf(m)
    enddo
 
    Sinf=kB*log(Zinf)
 
    write(20,*) "Infinite Temperature Limit"
    write(20,*) "        m  omega(m)   Erel(m)/eV       p(m)      p(m)/omega(m)"
    do m=1,Mm
       write(20,101) m, omega(m), Erel(m), pinf(m), pinf(m)/omega(m)
    enddo
 


    write(21,300)  "Infinite", Einf, " - ",Sinf
   300         format(a10,2x,f14.4,8x,a3,7x,e12.6)

    if (DATA_exists) then
       avedatainf(1:ncol)=0.0
       do m=1,Mm
          do col=1,ncol
             avedatainf(col)=avedatainf(col)+data(col,m)*pinf(m)
          enddo
       enddo
       write(22,301) adjustr("Infinite"), avedatainf(1:ncol) 
   301         format(a10,2x,10(f10.4,2x))
   endif

    if (SPECTRA_exists) then
       avespecinf(1:npoints)=0.0
       do point=1,npoints
          do m=1,Mm
             avespecinf(point)=avespecinf(point)+spec(point,m)*pinf(m)
          enddo
          if (avespecinf(point)/maxspec.lt.tolminspec) then
             avespecinf(point)=0.0
          endif
       enddo
   endif




    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !       Writing ave_spectra for all temperatures
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    if (SPECTRA_exists) then
       do point=1,npoints
          write(23,302) xspec(point),avespec(point,1:Ntt),avespecinf(point)
       enddo
   302         format(1(f10.3,2x),7(e12.6,2x))
    endif



    CLOSE(20)
    CLOSE(21)
    if (DATA_exists) CLOSE(22)
    if (SPECTRA_exists) CLOSE(23)


end

