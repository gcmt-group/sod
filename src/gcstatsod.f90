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

       PROGRAM gcstats 

       IMPLICIT NONE

       INTEGER,PARAMETER :: NCONFMAX=10000, NCOLMAX=10, NPOINTSMAX=800, NTEMPMAX=50, NOUTSODSMAX=30, NBISMAX=1000
       REAL (kind=8),PARAMETER :: kB=8.61734E-5, tolprob=1.0E-12, tolminspec=1.0E-6, eVA3toGPa=160.2176621
       INTEGER :: m, auxm, g, ncol, col, tt, Ntt, nsubs, nbis, npoints, point
       INTEGER :: mEmin, mEmax, nsubsEmin
       REAL (kind=8)  :: Emin, Emax, paux, maxspec
       INTEGER :: nsubsmin, nsubsmax
       REAL (kind=8),DIMENSION(NTEMPMAX) :: Z,E,F,S,T,Gpot
       REAL (kind=8) :: Zinf, Einf, Sinf
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NCONFMAX,NCOLMAX) :: data
       REAL (kind=8),DIMENSION(NCOLMAX,NTEMPMAX) :: avedata
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NCONFMAX,NPOINTSMAX) :: spec
       REAL (kind=8),DIMENSION(NPOINTSMAX,NTEMPMAX) :: avespec
       REAL (kind=8),DIMENSION(NCOLMAX) :: avedatainf
       REAL (kind=8),DIMENSION(NPOINTSMAX) :: xspec, avespecinf

       LOGICAL :: TEMPERATURES_exists, DATA_exists, INGC_exists, SPECTRA_exists, booldata, boolspec
       REAL (kind=8) :: xormuvalue, x, xeq, mu, nx, deltamu, oldmu 
       CHARACTER :: xormu*2
       CHARACTER :: charnsubs*1, filenameout*9, auxstring*15,filenameene*11
       CHARACTER,DIMENSION(0:NOUTSODSMAX) :: filenamedat*9, filenamespec*12
       INTEGER :: n, noutsods, nsubsread, npos
       INTEGER,DIMENSION(0:NOUTSODSMAX) :: Mm
       INTEGER,DIMENSION(0:NOUTSODSMAX,NCONFMAX) :: omega
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NCONFMAX)  :: ene, enemun, enemunrel 
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NCONFMAX,NTEMPMAX)  :: p
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NTEMPMAX)  :: pn
       CHARACTER(len=30) fmtemplist

       INTEGER :: omegasum, k
       REAL (kind=8)  :: naver, eneaver, epsilonA, epsilonB, epsilon, E0
       REAL (kind=8)  :: fx, fpx, A0, A1, A2, B0, B1, B2
       REAL (kind=8)  :: momentaA0, momentaA1, momentaA2, momentaB0, momentaB1, momentaB2
       REAL (kind=8)  :: a11, a12, a21, a22, bb1, bb2, alpha, beta, pbinomial
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX)  :: Qn, C, tau
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX,NCONFMAX)  :: pinf
       REAL (kind=8),PARAMETER :: tolmu=1.0E-10, tolq=1.0E-10
       
       REAL (kind=8) SUM, SUMSQ, valpolynom, qtest
       REAL (kind=8) valpolynomnew, valpolynomb, r1, r2, qa, qb, q, va, vb

       REAL (kind=8) lambda, v0, v1, bv, bm0, bm1, bb, vol, bm, voln, xn
       REAL (kind=8),DIMENSION(0:NOUTSODSMAX)  :: EVSC

!Input files

       INQUIRE(FILE="TEMPERATURES", EXIST=TEMPERATURES_exists)
       if (TEMPERATURES_exists) then
          OPEN (UNIT=10, FILE="TEMPERATURES")
       else
          write(*,*)
          write(*,*) "TEMPERATURES files not found - analysis will be done at 300 K and 1000 K only."
          write(*,*)
       endif

       INQUIRE(FILE="INGC", EXIST=INGC_exists)
       if (INGC_exists) then
          OPEN (UNIT=14, FILE="INGC")
       else
          write(*,*) "INGC file does not exist, but it is needed for grandcanonical analysis. Aborting execution."
          STOP 
       endif


!Output files
       OPEN (UNIT=20,FILE="probabilities.dat")
       OPEN (UNIT=21,FILE="thermodynamics.dat")


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
           CLOSE(10)
        else
           T(1)=300.0
           T(2)=1000.0
           Ntt=2
        endif



    write(*,*) "Performing grand-canonical analysis..."
    write(*,*) 
    write(*,*) 

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read the INGC file
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
    READ (14,*)
    READ (14,*) nsubsmin, nsubsmax
    READ (14,*)
    READ (14,*) xormu, xormuvalue
    READ (14,*)
    READ (14,*) lambda
    if (lambda > 0 ) then 
       READ (14,*)
       READ (14,*) v0, v1, bv
       READ (14,*)
       READ (14,*) bm0, bm1, bb
    endif
    CLOSE (14)

    if (xormu.eq."x ") then
        x=xormuvalue
    else
       mu=xormuvalue
    endif

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read the OUTSOD_xx files
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
    noutsods =  nsubsmax-nsubsmin+1

    do nsubs = nsubsmin, nsubsmax
       if (nsubs < 10 ) then 
           write (filenameout, "(a8,I1)") "OUTSOD_0",nsubs
           OPEN (UNIT=100+nsubs,FILE=trim(filenameout), status='old')
       else
           write (filenameout, "(a7,I2)") "OUTSOD_",nsubs
           OPEN (UNIT=100+nsubs,FILE=trim(filenameout), status='old')
       endif
    enddo 
    

    do nsubs = nsubsmin, nsubsmax
   
       read(100+nsubs,*) nsubsread, auxstring, auxstring, npos
       if (nsubsread.ne.nsubs) then
           write(*,*) 'ERROR: nsubsread.ne.nsubs'
           stop
       endif
   
       read(100+nsubs,*) Mm(nsubs)
       write(*,*) "Number of inequivalent configurations with ",nsubs," substitutions: ",Mm(nsubs)
       do m=1, Mm(nsubs)
         read(100+nsubs,*)  auxm, omega(nsubs,m)
       enddo
       CLOSE(100+nsubs)
       
    enddo
   
    !write(*,*) "nsubsmin,nsubsmax,npos=",nsubsmin,nsubsmax,npos
    !write(*,*) "nsubsmin/npos=",REAL(nsubsmin)/REAL(npos)
    !write(*,*) "nsubsmax/npos=",REAL(nsubsmax)/REAL(npos)
    if ((xormu.eq."x ").AND.((x.lt.REAL(nsubsmin)/REAL(npos)).or.(x.gt.REAL(nsubsmax)/REAL(npos)))) then
       write(*,*) "x is out of range, given the number of substitutions considered"
       write(*,*) "x should be between",REAL(nsubsmin)/REAL(npos),"and",REAL(nsubsmax)/REAL(npos)
       write(*,*) 
       STOP
    endif


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read ENERGIES_xx files 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
    do nsubs = nsubsmin, nsubsmax
       if (nsubs < 10 ) then
           write (filenameene, "(a10,I1)") "ENERGIES_0",nsubs
           OPEN (UNIT=200+nsubs,FILE=trim(filenameene), status='old')
       else
           write (filenameene, "(a9,I2)") "ENERGIES_",nsubs
           OPEN (UNIT=200+nsubs,FILE=trim(filenameene), status='old')
       endif
    enddo
   
    do nsubs = nsubsmin, nsubsmax
       do m=1,Mm(nsubs)
          read(200+nsubs,*) ene(nsubs,m)
       enddo
       CLOSE(200+nsubs)
    enddo

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read DATA_xx files 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   do nsubs = nsubsmin, nsubsmax
       if (nsubs < 10 ) then
           write (filenamedat(nsubs), "(a6,I1)") "DATA_0",nsubs
       else
           write (filenamedat(nsubs), "(a5,I2)") "DATA_",nsubs
       endif
    enddo

    DATA_exists=.TRUE. 
    do nsubs = nsubsmin, nsubsmax
       INQUIRE(FILE=filenamedat(nsubs), EXIST=booldata)
       if (booldata) then
          write(*,*) "DATA file found for", nsubs, "substitutions"
          OPEN (UNIT=300+nsubs,FILE=trim(filenamedat(nsubs)), status='old')
          read(300+nsubs,*) ncol
          do m=1,Mm(nsubs)
             read(300+nsubs,*) data(nsubs,m,1:ncol)
          enddo
          CLOSE (300+nsubs)
       endif
       DATA_exists=(booldata.and.DATA_exists)
    enddo

    if (DATA_exists) then
       write(*,*) 
       write(*,*) "Average data will be calculated."
       write(*,*) 
       OPEN (UNIT=22,FILE="ave_data.dat")
    else
       write(*,*) 
       write(*,*) "At least one DATA file not found. No average data will be calculated."
       write(*,*) 
    endif

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Read SPECTRA_xx files 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   do nsubs = nsubsmin, nsubsmax
       if (nsubs < 10 ) then
           write (filenamespec(nsubs), "(a9,I1)") "SPECTRA_0",nsubs
       else
           write (filenamespec(nsubs), "(a8,I2)") "SPECTRA_",nsubs
       endif
    enddo

    SPECTRA_exists=.TRUE.
    do nsubs = nsubsmin, nsubsmax
       INQUIRE(FILE=filenamespec(nsubs), EXIST=boolspec)
       if (boolspec) then
          write(*,*) "SPECTRA file found for", nsubs, "substitutions"
          open (UNIT=400+nsubs,FILE=trim(filenamespec(nsubs)), status='old')
          read(400+nsubs,*) npoints
          do m=1,Mm(nsubs) 
             read(400+nsubs,*) spec(nsubs,m,1:npoints)
          enddo
          CLOSE(400+nsubs)
       endif
       SPECTRA_exists=(boolspec.and.SPECTRA_exists)
    enddo

    if (SPECTRA_exists) then
       write(*,*) 
       write(*,*) "Average spectra will be calculated."
       write(*,*) 
       OPEN (UNIT=23,FILE="ave_spectra.dat")
       OPEN (UNIT=15, FILE="XSPEC")
       do point=1,npoints
          read(15,*)  xspec(point)
       enddo
       maxspec=MAXVAL(spec)
       CLOSE (15)
    else
       write(*,*) 
       write(*,*) "At least one SPECTRA file not found. No average spectra will be calculated."
       write(*,*) 
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
    !       This starts the loop over all temperature values
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    do tt=1,Ntt  ! loop over all temperature values
    
       write(*,*)
       write(*,*) "__________________________________________"
       write(*,*)
       write(*,*)             " T = ", T(tt), "K"
       write(*,*) "__________________________________________"
       write(*,*)
    


       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      If x is specified, get mu
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (xormu.eq."x ") then  ! if for xormu.eq."x "

!xx1

          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate omegasum
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          omegasum = 0
   
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                omegasum = omegasum + omega(nsubs,m)
             enddo
          enddo
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate naver
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          naver = 0.0 
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                naver = naver + omega(nsubs,m) * nsubs
             enddo
          enddo
          naver = naver / omegasum
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate eneaver
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          eneaver = 0.0
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                eneaver = eneaver + omega(nsubs,m) * ene(nsubs,m)
             enddo
          enddo
          eneaver = eneaver / omegasum
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate epsilon
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          epsilonA = 0.0
          epsilonB = 0.0
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                epsilonA = epsilonA + omega(nsubs,m) * (nsubs - naver) * (ene(nsubs,m) - eneaver)
                epsilonB = epsilonB + omega(nsubs,m) * (nsubs - naver) * (nsubs - naver)
             enddo
          enddo
          epsilon = epsilonA / epsilonB
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate E0
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          E0 = eneaver - epsilon * naver
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate EVSC (Energy of Volume-Stress Correction) 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          if (lambda == 0.0) then 
             write(*,*) "No Volume-Stress Correction applied."
             do nsubs = nsubsmin, nsubsmax
                EVSC(nsubs) = 0.0
             enddo
          else
             ! Calculate Volume at this x
             vol = v0  * (1-x) + v1  * x + bv * x * (1-x)
             bm  = (bm0 * (1-x) + bm1 * x + bb * x * (1-x)) / eVA3toGPa
             do nsubs = nsubsmin, nsubsmax
                xn = REAL(nsubs)/REAL(npos)
                voln = v0  * (1-xn) + v1  * xn + bv * xn * (1- xn)
                EVSC(nsubs) = (lambda / 2.0) * bm * vol * (1 - (voln/vol))**2
                write(*,*) "Energy of Volume-Stress Correction: nsubs=",nsubs,"voln=",voln,"EVSC(nsubs)=",EVSC(nsubs) 
             enddo
          endif


          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate Qn 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do nsubs = nsubsmin, nsubsmax
             Qn(nsubs) = 0.0
             do m=1,Mm(nsubs)
                Qn(nsubs) = Qn(nsubs) + REAL(omega(nsubs,m))*exp(-(ene(nsubs,m)+EVSC(nsubs)-E0-nsubs*epsilon)/(kB*T(tt)))
             enddo
          enddo
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate polynomial coefficients Cn
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do nsubs = nsubsmin, nsubsmax
             C(nsubs) = (real(nsubs)/real(npos)-x) * Qn(nsubs)
          enddo
   
   
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !      Calculate q with the bisection method
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write(*,*) "Using bisection method to find chemical potential corresponding to specified composition"
          q = x / (1 - x)
          mu = epsilon + kB * T(tt) * LOG(q)
          write(*,*) "Initial chemical potential:  mu =", mu, " eV"

          valpolynom = 0.0
          do nsubs = nsubsmin, nsubsmax
             valpolynom = valpolynom + C(nsubs) * (q**nsubs)
          enddo

          valpolynomnew=valpolynom
          do while (valpolynom * valpolynomnew > 0 )

             call random_seed()
             call random_number(r1)
             call random_seed()
             call random_number(r2)
             do while (r1 < tolq) 
                 call random_seed()
                 call random_number(r1)
             enddo
             do while (r2 < tolq) 
                 call random_seed()
                 call random_number(r2)
             enddo

             qtest = q * r1/r2
             valpolynomnew = 0.0
             do nsubs = nsubsmin, nsubsmax
                valpolynomnew = valpolynomnew +  C(nsubs) * (qtest**nsubs)
             enddo

          enddo

          if (q < qtest) then
             qa = q 
             qb = qtest
             va = valpolynom
             vb = valpolynomnew
          else 
             qa = qtest
             qb = q
             va = valpolynomnew
             vb = valpolynom
          endif

          do nbis = 1, NBISMAX 

             q = (qa + qb)/2
             oldmu = mu
             mu = epsilon + kB * T(tt) * LOG(q)
             valpolynomnew = 0.0
             do nsubs = nsubsmin, nsubsmax
                valpolynomnew = valpolynomnew + C(nsubs) * (q**nsubs)
             enddo

             if (valpolynomnew * va > 0) then
                qa = q
                va = valpolynomnew
             else
                qb = q
                vb = valpolynomnew
             endif

             if (( ABS(qa-qb) < tolq ).AND.(ABS(mu-oldmu)< tolmu))  then
                write(*,*) "Convergence achieved:  delta(mu) = ", mu-oldmu, " eV"
                EXIT
             endif

          enddo


          mu = epsilon + kB * T(tt) * log(q)
          write(*,*) "Converged chemical potential: mu = ",mu, " eV"

          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          Get minimum value of the grand potential to calculate enemunrel
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                enemun(nsubs,m) = ene(nsubs,m)+EVSC(nsubs)-nsubs*mu
             enddo
          enddo
         
          Emin=enemun(nsubsmin,1)
          mEmin=1
          nsubsEmin=nsubsmin
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                if (Emin > enemun(nsubs,m)) then
                    Emin=enemun(nsubs,m)
                    mEmin=m
                    nsubsEmin=nsubs
                endif
             enddo
          enddo
         
          do nsubs = nsubsmin, nsubsmax
              do m=1,Mm(nsubs)
                 enemunrel(nsubs,m) = enemun(nsubs,m) - Emin
              enddo
          enddo
         
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          Calculate the partition function from mu
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          Z(tt)=0.0
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                  Z(tt)=Z(tt)+omega(nsubs,m)*exp(-enemunrel(nsubs,m)/(kB*T(tt)))
             enddo 
          enddo
          
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          Calculate probabilities and x (equilibrium composition) from mu 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          nx = 0
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                p(nsubs,m,tt)=omega(nsubs,m)*exp(-enemunrel(nsubs,m)/(kB*T(tt)))/Z(tt)
                nx = nx + nsubs*p(nsubs,m,tt)
             enddo
          enddo
          
          xeq= nx/npos
          write(*,*) "Composition calculated for this chemical potential: x = ", xeq

!xx2
       else !(xormu = 'mu')
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          Get minimum value of the grand potential to calculate enemunrel
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                enemun(nsubs,m) = ene(nsubs,m)-nsubs*mu
             enddo
          enddo

          Emin=enemun(nsubsmin,1)
          mEmin=1
          nsubsEmin=nsubsmin
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                if (Emin > enemun(nsubs,m)) then
                    Emin=enemun(nsubs,m)
                    mEmin=m
                    nsubsEmin=nsubs
                endif
             enddo
          enddo

          do nsubs = nsubsmin, nsubsmax
              do m=1,Mm(nsubs)
                 enemunrel(nsubs,m) = enemun(nsubs,m) - Emin
              enddo
          enddo
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          Calculate the partition function from mu
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          Z(tt)=0.0
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                  Z(tt)=Z(tt)+omega(nsubs,m)*exp(-enemunrel(nsubs,m)/(kB*T(tt)))
             enddo
          enddo

          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !          Calculate probabilities and x (equilibrium composition) from mu 
          !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          nx = 0
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                p(nsubs,m,tt)=omega(nsubs,m)*exp(-enemunrel(nsubs,m)/(kB*T(tt)))/Z(tt)
                nx = nx + nsubs*p(nsubs,m,tt)
             enddo
          enddo

          xeq= nx/npos
          write(*,*) "Composition at this temperature, for the given chemical potential: x = ", xeq
          x = xeq

       endif ! End if for xormu.eq."x "
!xx3

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !          Writing  probabilities.txt file
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(20,*) 
       write(20,*) "____________________________________________________________________________________________"
       write(20,*) "Temperature: T = ",T(tt), "K"
       write(20,*) "Chemical potential: mu = ",mu, "eV"
       write(20,*) "Composition: x = ", xeq
       write(20,*) 
    
       nx = 0
       write(20,*) " n      m     omega(n,m)       E(n,m)          E-n*mu         prob(n,m,T)   prob(n,m,T)/omega    "
       do nsubs = nsubsmin, nsubsmax
          pn(nsubs,tt) = 0.0
          do m=1,Mm(nsubs)
             write(20,101) nsubs, m, omega(nsubs,m), ene(nsubs,m),  enemun(nsubs,m),  p(nsubs,m,tt), p(nsubs,m,tt)/omega(nsubs,m)
101          format(      i3,1x, i6,    1x,i8,5x,   2(4x,f14.5)                       2(4x,e12.6))
             pn(nsubs,tt) = pn(nsubs,tt) + p(nsubs,m,tt)
          enddo
       enddo

       write(20,*) "------------------------------"
       write(20,*) " n      Cumulative prob(n,T) "
       do nsubs = nsubsmin, nsubsmax
             write(20,102) nsubs,pn(nsubs,tt)
102          format(      i6,    5x,   f14.5 )
       enddo
       write(20,*) "------------------------------"

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      Calculate the energy (E), free energy (F), entropy (S) 
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
       E(tt)=0.0
       S(tt)=0.0
      
       do nsubs = nsubsmin, nsubsmax
          do m=1,Mm(nsubs)
             E(tt)=E(tt)+ene(nsubs,m)*p(nsubs,m,tt)
             S(tt)=S(tt)-kB*p(nsubs,m,tt)*log(p(nsubs,m,tt)/omega(nsubs,m))
          enddo
       enddo
      
       F(tt) = E(tt) - T(tt)*S(tt)
      
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      Calculate the average data
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
       if (DATA_exists) then
        
           avedata(1:ncol,tt)=0.0
        
           do col=1,ncol
               do nsubs = nsubsmin, nsubsmax
                  do m=1,Mm(nsubs)
                     avedata(col,tt)=avedata(col,tt)+data(nsubs,m,col)*p(nsubs,m,tt)
                  enddo
               enddo
           enddo
        
        endif

       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !      Calculate the average spectra
       !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       if (SPECTRA_exists) then

           avespec(1:npoints,tt)=0.0

           do point=1,npoints
              do nsubs = nsubsmin, nsubsmax
                 do m=1,Mm(nsubs)
                    avespec(point,tt)=avespec(point,tt)+spec(nsubs,m,point)*p(nsubs,m,tt)
                 enddo
              enddo
              if (avespec(point,tt)/maxspec.lt.tolminspec) then
                 avespec(point,tt)=0.0
              endif
           enddo

        endif

   
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
    !      Calculate the limit of full disorder
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Calculate residual momenta
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    A0 = momentaA0(nsubsmax,npos,x)
    A1 = momentaA1(nsubsmax,npos,x)
    A2 = momentaA2(nsubsmax,npos,x)
    B0 = momentaB0(nsubsmin,npos,x)
    B1 = momentaB1(nsubsmin,npos,x)
    B2 = momentaB2(nsubsmin,npos,x)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Calculate alpha and beta
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    a11 = 1.0 - B0 - A0
    a12 = x * npos - B1 - A1
    a21 = a12
    a22 = ( x * npos * ( 1.0 + x * (npos -1) ) ) - B2 - A2
    bb1 = 1.0
    bb2 = x * npos

    alpha = ( (bb1*a22) - (a12*bb2) ) / ( (a11*a22) - (a12*a21) )
    beta  = ( (a11*bb2) - (a21*bb1) ) / ( (a11*a22) - (a12*a21) )

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Calculate tau
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    do nsubs = nsubsmin, nsubsmax
        tau(nsubs) = alpha + beta*nsubs
    enddo

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !      Calculate renormalised probabilities and equilibrium concentration
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    nx=0.0
    do nsubs = nsubsmin, nsubsmax
       do m=1,Mm(nsubs)
          pinf(nsubs,m) = tau(nsubs) * omega(nsubs,m) * (x**nsubs) * ((1-x)**(npos-nsubs))
          nx = nx + nsubs*pinf(nsubs,m)
       enddo
    enddo    
 
 
    xeq= nx/npos

    write(*,*)
    write(*,*) "___________________________________________________"
    write(*,*)
    write(*,*)             " T = infinity (ideal disorder limit)"
    write(*,*) "___________________________________________________"
    write(*,*)
    write(*,*) "Composition from probabilities:"
    write(*,*) "x = ", xeq


    if (xormu.eq."mu") then
       x = xeq
    endif



    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !          Write probabilities in the limit of infinite temperature
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    write(20,*)
    write(20,*) "____________________________________________________________________________________________"
    write(20,*) "Ideal disorder limit"
    write(20,*) "Composition: x = ", xeq
    write(20,*) " n      m     omega(n,m)     E(n,m)                           prob(n,m,T)   prob(n,m,T)/omega    "
    do nsubs = nsubsmin, nsubsmax
       do m=1,Mm(nsubs)
            write(20,103) nsubs, m,    omega(nsubs,m), ene(nsubs,m),   pinf(nsubs,m), pinf(nsubs,m)/omega(nsubs,m)
103         format(       i3, 1x,i6,1x,i8,5x,          4x,f14.5,   16x,2(4x,e12.6))
       enddo
    enddo


    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !          Calculate and write results in the limit of infinite temperature
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    Einf=0.0
    do nsubs = nsubsmin, nsubsmax
       do m=1,Mm(nsubs)
          Einf=Einf+ene(nsubs,m)*pinf(nsubs,m)
       enddo
    enddo

    Sinf= -kB * ( x*log(x)+(1-x)*log(1-x) ) * npos

    ! Write thermodynamics.dat
    write(21,300)  "Infinite", Einf, " - ",Sinf
   300         format(a10,2x,f14.4,8x,a3,7x,e12.6)



    ! Calculate and write avedatainf 
    if (DATA_exists) then
        avedatainf(1:ncol)=0.0
        do col=1,ncol
            do nsubs = nsubsmin, nsubsmax
                do m=1,Mm(nsubs)
                    avedatainf(col)=avedatainf(col)+data(nsubs,m,col)*pinf(nsubs,m)
                enddo
            enddo
        enddo
       write(22,301) adjustr("Infinite"), avedatainf(1:ncol)
   301         format(a10,2x,10(f10.4,2x))
     endif



    ! Calculate avespecinf 
    if (SPECTRA_exists) then
       avespecinf(1:npoints)=0.0
       do point=1,npoints
          do nsubs = nsubsmin, nsubsmax
             do m=1,Mm(nsubs)
                avespecinf(point)=avespecinf(point)+spec(nsubs,m,point)*pinf(nsubs,m)
             enddo
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
302       format(1(f10.3,2x),8(e12.6,2x))
       enddo
    endif

   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    CLOSE(20)
    CLOSE(21)
    if (DATA_exists) CLOSE(22)
    if (SPECTRA_exists) CLOSE(23)

    write(*,*)
    write(*,*)
    write(*,*) "---------------------------------------"
    write(*,*) "Grand-canonical analysis completed."
    write(*,*) "---------------------------------------"
    write(*,*)


end

