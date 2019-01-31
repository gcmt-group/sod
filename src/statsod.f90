!*******************************************************************************
!    Copyright (c) 2019 Ricardo Grau-Crespo, Said Hamad
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

       INTEGER,PARAMETER :: NCONFMAX=100000, NCOLMAX=20, NTEMPMAX=100
       REAL*8,PARAMETER :: kB=8.61734E-5
       INTEGER :: m, auxm, Mm, g, ncol, col, tt, Ntt, nsubs
       INTEGER, DIMENSION(NTEMPMAX) :: mpmax, mpmin
       REAL*8,DIMENSION(NTEMPMAX) :: Z,E,F,S,Edmin,Edmax,T
       REAL*8,DIMENSION(NCONFMAX) :: ene, Sd
       REAL*8,DIMENSION(NCONFMAX,NTEMPMAX) :: p,Ed,Edrel
       INTEGER,DIMENSION(NCONFMAX) :: conf,omega
       REAL*8,DIMENSION(NCOLMAX,NCONFMAX) :: data 
       REAL*8,DIMENSION(NCOLMAX,NTEMPMAX) :: avedata 
       LOGICAL :: TEMPERATURES_exists, DATA_exists



!Input files
       OPEN (UNIT=11,FILE="OUTSOD")
       OPEN (UNIT=12,FILE="ENERGIES")
       OPEN (UNIT=13,FILE="DATA")

       INQUIRE(FILE="TEMPERATURES", EXIST=TEMPERATURES_exists)
       if (TEMPERATURES_exists) then
          OPEN (UNIT=10, FILE="TEMPERATURES")
       endif


!Output files
       OPEN (UNIT=20,FILE="probabilities.dat")
       OPEN (UNIT=21,FILE="statistics.dat")


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
        else
           T(1)=1.0
           T(2)=300.0
           T(3)=1000.0
           T(4)=100000000000.0
           Ntt=4
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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Read ENERGIES 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do m=1,Mm
	read(12,*) ene(m)
	enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Read the data file, DATA 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        read(13,*) ncol
	do m=1,Mm
	read(13,*) data(1:ncol,m) 
	enddo


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Start writing STATISTICS file
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        write(21,*) "      T        H           G          S             average values of DATA" 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       This starts the loop over all temperature values
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do tt=1,Ntt


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Get reference values for Ed (its minimum value at its temperature), 
!          in order to get better accuracy for the exponential calculations
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


           Sd(1)=kB*log(REAL(omega(1)))
	   Ed(1,tt)=ene(1)-T(tt)*Sd(1)
 	   Edmin(tt)=Ed(1,tt)
 	   Edmax(tt)=Ed(1,tt)
           mpmax(tt)=1
           mpmin(tt)=1
           do m=2,Mm
              Sd(m)=kB*log(REAL(omega(m)))
	      Ed(m,tt)=ene(m)-T(tt)*Sd(m)
              if (Edmin(tt) > Ed(m,tt)) then
                  Edmin(tt)=Ed(m,tt)
	          mpmax(tt)=m
              endif 
              if (Edmax(tt) < Ed(m,tt)) then
                  Edmax(tt)=Ed(m,tt)
	          mpmin(tt)=m
              endif 
	   enddo			
           
	   write(20,110) T(tt)
110        format(" T = ", f6.1, " K")
	   write(20,*) "Configuration with maximum probability: ",mpmax(tt)
	   write(20,*) "Configuration with minimum probability: ",mpmin(tt)
	   write(20,*) 


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Calculate the partition function and probabilities
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  	   Z(tt)=0.0
	   do m=1,Mm
	      Edrel(m,tt)=Ed(m,tt)-Edmin(tt)
	      Z(tt)=Z(tt)+exp(-Edrel(m,tt)/(kB*T(tt)))
	   enddo										  

	   do m=1,Mm
	      p(m,tt)=exp(-Edrel(m,tt)/(kB*T(tt)))/Z(tt)	
	   enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the energy (E) and Helmholtz energy (F)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   E(tt)=0.0
	   S(tt)=0.0
           avedata(1:ncol,tt)=0.0
	 
	   do m=1,Mm
	      E(tt)=E(tt)+ene(m)*p(m,tt)
	      do col=1,ncol
	      avedata(col,tt)=avedata(col,tt)+data(col,m)*p(m,tt)
	      enddo
	   enddo

	   F(tt)=Edmin(tt)-kB*T(tt)*log(Z(tt))

	   S(tt)=(E(tt)-F(tt))/T(tt)
        

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Write PROBABILITIES for temperature T 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   
           write(20,*) "        m  omega(m)         Sd(m)      ene(m)       Ed(m)        p(m)"
	   do m=1,Mm
	      write(20,100) m,omega(m),Sd(m),ene(m),Ed(m,tt),p(m,tt)
  100         format(i10,2x,i8,2x,e12.6,3(2x,f10.4))
	   enddo

	   write(20,*) 
	   write(20,*)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!          Write STATISTICS for temperature T 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           write(21,200)  T(tt), E(tt), F(tt), S(tt), avedata(1:ncol,tt) 
  200      format(f10.1,2x,2(f10.4,2x),e12.6,2x,9(f10.4,2x))

	enddo

	end

