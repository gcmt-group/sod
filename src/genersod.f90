
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

    PROGRAM genersod 
       IMPLICIT NONE

       INTEGER,PARAMETER :: NSPMAX=10 , NATMAX=10000, NOPMAX=5000, NCELLMAX=1000
       INTEGER,PARAMETER :: NLINEAMAX=200 
       REAL,PARAMETER :: tol0=0.0001, kB=8.61734E-5
      
       INTEGER :: i,j,k,l,xx,yy,suma,t,ina,inb,inc,elei, factorial, sub, m
       INTEGER :: op1, nop1, op, nop, opsc, nopsc,ifound,op1new,nop1new,aux,opc,nopc
       INTEGER :: sp, nsp, spmap, nspmap, cumnatsp, cumnatspmap, sptarget, sptargetmap,ssp
       INTEGER :: at0, nat0, at1, nat1, nat1r, at, nat , atmap, natmap, at1r, at1i,attmp,att,mapno,FILER, MAPPER
       INTEGER :: na,nb,nc,nsubs,nsubsmap,atini,atfin,atinimap,atfinmap
       INTEGER :: pos,npos,count,ntc,ntcmax,nic,equivcount,iequiv,indcount
       LOGICAL :: found,foundnoind,mores
       INTEGER,DIMENSION(:), ALLOCATABLE:: newconf
       INTEGER,DIMENSION(2)   :: newshell, newshellmap
       INTEGER,DIMENSION(NOPMAX)   :: op1good
       INTEGER,DIMENSION(:), ALLOCATABLE:: degen
       INTEGER,DIMENSION(NSPMAX) :: natsp0, natsp1, natsp, natspmap, snatsp, ishell, ishellmap
       INTEGER,DIMENSION(NATMAX) :: spat0, spat1, spat, spatmap, spat1r
       REAL,DIMENSION(NATMAX,3) :: coords0, coords1, coords, coords1r, result, coordsmap
       REAL,DIMENSION(3) :: coordstemp
       REAL,DIMENSION(3,3) :: cellvector
       REAL,DIMENSION(NOPMAX,3,3) :: mgroup1, mgroup, mgroup1new
       REAL,DIMENSION(NOPMAX,3)   :: vgroup1, vgroup,vgroup1sca,vgroup1scb,vgroup1scc,vgroup1new
       REAL,DIMENSION(NCELLMAX,3)   :: vt
       REAL,DIMENSION(3) :: x, vcorr
       INTEGER,DIMENSION(NOPMAX,NATMAX)   :: fulleqmatrix,eqmatrixtarget
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: conf,indconf,indconfmap,equivconf
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: mapping
       INTEGER,DIMENSION(NATMAX)   :: as,test
       REAL :: a1,b1,c1,alpha,beta,gamma,a,b,c,cc, amap, bmap, cmap, alphamap, betamap, gammamap
       REAL :: prod, xs, ientropyguess, ientropy, perc
       CHARACTER,DIMENSION(NSPMAX) :: symbol*3, symbolmap*3, ssymbol*3
       CHARACTER,DIMENSION(2) :: newsymbol*3, newsymbolmap*3
       CHARACTER :: linea*85, title*15, tmptitle*15, runtitle*40, trashtext*20, symboltrash*3
       CHARACTER :: sindcount*5,snumber*4

! Input files

       OPEN (UNIT= 9,FILE="INSOD")
       OPEN (UNIT=30,FILE="OUTSOD")
       OPEN (UNIT=31,FILE="SUPERCELL")

! Output files

       OPEN (UNIT=43,FILE="filer")

!
! DEFINITION OF VARIABLES:
!
! sp                  Index for the species 
! nsp                 Total number of species
! sptarget            Number of the species to be substituted
! at0,at1,at          Indexes for the atoms in the assymetric unit, unit cell and supercell
! nat0,nat1,nat       Total numbers of atoms in the assymetric unit, unit cell and supercell
! at1r,nat1r	        Idem for the atoms in the redundant cell (with repeated positions)
! atini,atfin	        Initial and final atom indexes of the species to be substituted
! pos                 Index for atomic positions of the target species
! npos		            Number of atoms of the target species 
! conf		            List of all configurations (each configuration is a list of the substituted positions)
! count		            Index for the configurations (conf)     
! ntc                 Total number of configurations in conf (count=1,ntc)
! indconf             List of independent configurations 
! indcount            Index for the independent configurations
! nic                 Total number of independent configurations in indconf (indcount=1,nic)
! equivconf           Temporary list containing the equivalent configurations at every step of the algorithm
! equivcount          Index for the equivalent configurations (equivconf)     
! tol0 		            General tolerance
! tol1 		            Tolerance used for correcting the x-FLOOR(x) function
!
!
!

       WRITE (*,*) "**************************************************************************** " 
       WRITE (*,*) "         SOD (Site Occupancy Disorder) version 0.51  " 
       WRITE (*,*) " " 
       WRITE (*,*) "         Authors: R. Grau-Crespo and S. Hamad                                   " 
       WRITE (*,*) " " 
       WRITE (*,*) "         Contact:  <r.grau-crespo@reading.ac.uk> " 
       WRITE (*,*) "**************************************************************************** " 
       WRITE (*,*) " " 
       WRITE (*,*) " " 
       WRITE (*,*) " " 
       WRITE (*,*) "Reading INSOD, SUPERCELL, and OUTSOD to generate calculation input files" 
       WRITE (*,*) " " 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the INSOD file 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       READ (9,*)

       READ (9,*) runtitle
       READ (9,*)
       READ (9,*)
       READ (9,*) a1,b1,c1,alpha,beta,gamma
       READ (9,*)
       READ (9,*)
       READ (9,*) nsp
       READ (9,*)
       READ (9,*)
       READ (9,*)(symbol(sp),sp=1,nsp)
       READ (9,*)
       READ (9,*)
       READ (9,*)(natsp0(sp),sp=1,nsp)
       READ (9,*)
       READ (9,*)

       nat0=0
       do sp=1,nsp
         nat0 = nat0 + natsp0(sp)
       enddo

       do at0=1,nat0
         READ (9,*)(coords0(at0,i),i=1,3)
       enddo
       READ (9,*)
       READ (9,*)
       READ (9,*) na,nb,nc
       READ (9,*)
       READ (9,*)
       READ (9,*) sptarget
       READ (9,*)
       READ (9,*)
       READ (9,*) nsubs
       If (nsubs ==0) then
        Write(*,*) "Illegal number of substitutions"
        Stop
       Endif
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*) (newsymbol(i),i=1,2)
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*)
       READ (9,*) FILER, MAPPER
       READ (9,*)

       if ((FILER>0).AND.(FILER<10)) then
       	READ (9,*)
       	READ (9,*)
       	READ (9,*) (ishell(sp),sp=1,nsp)
       	READ (9,*)
       	READ (9,*) (newshell(i),i=1,2)
       endif




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the OUTSOD file 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	ntcmax=10000

        READ(30,*) nsubs, trashtext, trashtext, npos
        READ(30,*) nic 

!!!!!!!!Allocating array sizes

        ALLOCATE(degen(1:nic))
        ALLOCATE(newconf(1:nsubs))
        ALLOCATE(conf(1:nic,1:nsubs))
        ALLOCATE(indconf(1:nic,1:nsubs))

        do indcount=1,nic

                READ(30,*) m,degen(indcount),indconf(indcount,1:nsubs)
                if (m .ne. indcount) then 
                        WRITE (*, *) "Error in configuration numbering in OUTSOD. Aborting..."
                        STOP
                Endif
        enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the SUPERCELL file 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        READ(31,*) a,b,c,alpha,beta,gamma
        READ(31,*) natsp(1:nsp)

        nat=sum(natsp)
        do at=1,nat
                 READ(31,*) symboltrash,coords(at,1),coords(at,2),coords(at,3)
        enddo


!cccccccccccccccccccccccccccccccccccc
! Generating spat array
!cccccccccccccccccccccccccccccccccccc

       do at=1, natsp(1)
          spat(at)=1
       enddo 
       cumnatsp=natsp(1)
       do sp=2, nsp
          do at=cumnatsp+1, cumnatsp+natsp(sp)
             spat(at)=sp
          enddo
          cumnatsp=cumnatsp+natsp(sp) 
       enddo



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Calculate the initial and final nat of the target species
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if (sptarget.eq.1) then 
                atini=1
        else
                atini=1
                do sp=1,sptarget-1
                        atini=atini+natsp(sp)
                enddo
        endif      

        atfin=atini+natsp(sptarget)-1



!!!!!!!!Allocating array sizes

!        ALLOCATE(degen(1:nic))
!        ALLOCATE(newconf(1:nsubs))
!        ALLOCATE(conf(1:nic,1:nsubs))
!        ALLOCATE(indconf(1:nic,1:nsubs))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!    GENERATE INPUT FILES !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! Write a temporary file with FILER, to be read later by the script
	WRITE(43,*) FILER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SELECT CASE(FILER)


        CASE(0)
                WRITE(*,*) "Calculation files not created.&
                           & Change FILER value in INSOD file if you want to create calculation files."
        WRITE(*,*) ""
              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!	GULP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CASE(1)

        WRITE (*,*) " " 
        WRITE (*,*) "Creating input files for GULP in ./CALCS... " 
        WRITE (*,*) " " 

        OPEN (UNIT=41,FILE="top.gulp")
        OPEN (UNIT=42,FILE="bottom.gulp")

        do indcount=1,nic
	   newconf(1:nsubs)=indconf(indcount,1:nsubs)

	do l=1,NLINEAMAX
	        READ (41,332,end=99) linea
 332   		format(a85)
		WRITE(indcount+100000,332) linea
	enddo
  99    continue
	rewind(41)

        WRITE(indcount+100000,211) a,b,c,alpha,beta,gamma
	WRITE(indcount+100000,*) "frac"
 211    format(6(f10.4,2x))
        do at=1,nat
		sp=spat(at)
		if (sp.ne.sptarget) then
                    WRITE(indcount+100000,331) symbol(sp),"core",coords(at,1),coords(at,2),coords(at,3)
	            if (ishell(sp).eq.1) WRITE(indcount+100000,331) symbol(sp),"shel",coords(at,1),coords(at,2),coords(at,3)
		else
		    att=at-atini+1
                    call member(nsubs,newconf,att,ifound)
		    if (ifound==1) then
		         WRITE (indcount+100000,331) newsymbol(1),"core",coords(at,1),coords(at,2),coords(at,3)
	                 if (newshell(1).eq.1) WRITE(indcount+100000,331) newsymbol(1),"shel",coords(at,1),coords(at,2),coords(at,3)
                    else
		         WRITE (indcount+100000,331) newsymbol(2),"core",coords(at,1),coords(at,2),coords(at,3)
	                 if (newshell(2).eq.1) WRITE(indcount+100000,331) newsymbol(2),"shel",coords(at,1),coords(at,2),coords(at,3)
		    endif
		endif
 331    format(a3,2x,a4,2x,3(f11.7,2x))
        enddo


	do l=1,NLINEAMAX
		READ (42,332,end=199) linea
		WRITE(indcount+100000,333) linea
 333		format(a85)
	enddo
 199    continue
	rewind(42)

!        if (indcount<10) then
!        WRITE(UNIT=snumber, FMT='(I1)') indcount
!        WRITE (indcount+100000,334) "output xtl c0000",snumber
!        WRITE (indcount+100000,334) "output arc c0000",snumber
!        endif
!        if ((indcount>9).AND.(indcount<100)) then
!        WRITE(UNIT=snumber, FMT='(I2)') indcount
!        WRITE (indcount+100000,334) "output xtl  c000",snumber
!        WRITE (indcount+100000,334) "output arc  c000",snumber
!        endif
!        if ((indcount>99).AND.(indcount<1000)) then
!        WRITE(UNIT=snumber, FMT='(I3)') indcount
!        WRITE (indcount+100000,334) "output xtl   c00",snumber
!        WRITE (indcount+100000,334) "output arc   c00",snumber
!        endif
!        if ((indcount>999).AND.(indcount<10000)) then
!        WRITE(UNIT=snumber, FMT='(I4)') indcount
!        WRITE (indcount+100000,334) "output xtl    c0",snumber
!        WRITE (indcount+100000,334) "output arc    c0",snumber
!        endif
!        if ((indcount>9999).AND.(indcount<100000)) then
!        WRITE(UNIT=snumber, FMT='(I4)') indcount
!        WRITE (indcount+100000,334) "output xtl    c",snumber
!        WRITE (indcount+100000,334) "output arc    c",snumber
!        endif
!
!        if (indcount>99999) then
!	WRITE (indcount+100000,*) "Error, too many configurations (>99999)! Calculation files not written!"
!	endif
! 334    format(a15,a4)

        CLOSE (UNIT=indcount+100000)

        enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! METADISE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CASE(2)

        WRITE (*,*) " " 
        WRITE (*,*) "Creating input files for METADISE in ./CALCS... " 
        WRITE (*,*) " " 

        OPEN (UNIT=51,FILE="top.metadise")
        OPEN (UNIT=52,FILE="bottom.metadise")

        do indcount=1,nic
           newconf(1:nsubs)=indconf(indcount,1:nsubs)

        do l=1,NLINEAMAX
           READ (51,332,end=299) linea
           WRITE(indcount+100000,332) linea
        enddo
 299    continue
        rewind(51)

        WRITE(indcount+100000,212) "CELL ",a,b,c,alpha,beta,gamma
 212    format(a5,6(f10.4,2x))
        WRITE(indcount+100000,*) "SPACE FULL P1 1 1"
        WRITE(indcount+100000,*) "FRAC"
        do at=1,nat
                sp=spat(at)
                if (sp.ne.sptarget) then
                    WRITE(indcount+100000,331) symbol(sp),"CORE",coords(at,1),coords(at,2),coords(at,3)
                    if (ishell(sp).eq.1) WRITE(indcount+100000,331) symbol(sp),"SHEL",coords(at,1),coords(at,2),coords(at,3)
                else
                    att=at-atini+1
                    call member(nsubs,newconf,att,ifound)
                    if (ifound==1) then
                         WRITE (indcount+100000,331) newsymbol(1),"CORE",coords(at,1),coords(at,2),coords(at,3)
                         if (newshell(1).eq.1) WRITE(indcount+100000,331) newsymbol(1),"SHEL",coords(at,1),coords(at,2),coords(at,3)
                    else
                         WRITE (indcount+100000,331) newsymbol(2),"CORE",coords(at,1),coords(at,2),coords(at,3)
                         if (newshell(2).eq.1) WRITE(indcount+100000,331) newsymbol(2),"SHEL",coords(at,1),coords(at,2),coords(at,3)
                    endif
                endif
        enddo
        do l=1,NLINEAMAX
        READ (52,332,end=399) linea
        WRITE(indcount+100000,333) linea
        enddo
 399   continue
        rewind(52)

        CLOSE (UNIT=indcount+100000)

        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!	VASP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CASE(11)

        WRITE (*,*) " " 
        WRITE (*,*) "Creating input files for VASP in ./CALCS... " 
        WRITE (*,*) " " 


        title='vasp'
!        do sp=1,nsp
!           if (sp.ne.sptarget) then
!              tmptitle=title//symbol(sp)
!              title=tmptitle
!           else
!              tmptitle=title//newsymbol(1)
!              title=tmptitle//newsymbol(2)
!           endif
!        enddo
!
!	title=symbol(1)//symbol(2)

        do indcount=1,nic
        newconf(1:nsubs)=indconf(indcount,1:nsubs)

        WRITE(indcount+100000,*) title
        WRITE(indcount+100000,*) '1.00000000' 
        call cell(cellvector, a, b,c,alpha,beta,gamma)
        WRITE(indcount+100000,335) cellvector(1,1), cellvector(2,1), cellvector(3,1) 
        WRITE(indcount+100000,335) cellvector(1,2), cellvector(2,2), cellvector(3,2) 
        WRITE(indcount+100000,335) cellvector(1,3), cellvector(2,3), cellvector(3,3) 
 335    format(3(f10.6,2x))

        do ssp=1,nsp
           if (ssp<sptarget) then
              snatsp(ssp)=natsp(ssp)
              ssymbol(ssp)=symbol(ssp)
           endif
           if (ssp==sptarget) then
              snatsp(ssp)=nsubs
              snatsp(ssp+1)=npos-nsubs
              ssymbol(ssp)=newsymbol(1)
              ssymbol(ssp+1)=newsymbol(2)
           endif
           if (ssp>sptarget) then
              snatsp(ssp+1)=natsp(ssp)
              ssymbol(ssp+1)=symbol(ssp)
           endif
        enddo


        WRITE(indcount+100000,*) ssymbol(1:nsp+1) 
        WRITE(indcount+100000,336) snatsp(1:nsp+1)
 336    format(10(i4,1x))
        WRITE(indcount+100000,*) 'Direct'
        do at=1,atini-1
           sp=spat(at)
           WRITE(indcount+100000,335) coords(at,1),coords(at,2),coords(at,3)
	enddo
        do at=atini,atfin
           sp=spat(at)
           att=at-atini+1
           call member(nsubs,newconf,att,ifound)
           if (ifound==1) then
               WRITE (indcount+100000,335) coords(at,1),coords(at,2),coords(at,3)
           endif
	enddo
        do at=atini,atfin
           sp=spat(at)
           att=at-atini+1
           call member(nsubs,newconf,att,ifound)
           if (ifound==0) then
               WRITE (indcount+100000,335) coords(at,1),coords(at,2),coords(at,3)
           endif
        enddo
        do at=atfin+1,nat
           sp=spat(at)
           WRITE(indcount+100000,335) coords(at,1),coords(at,2),coords(at,3)
	enddo

        CLOSE (indcount+100000)

        enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!	CASTEP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CASE(12)

        WRITE (*,*) " " 
        WRITE (*,*) "Creating input files for CASTEP in ./CALCS... " 
        WRITE (*,*) " " 

        OPEN (UNIT=62,FILE="bottom.castep")

        do indcount=1,nic
   	   newconf(1:nsubs)=indconf(indcount,1:nsubs)
   
           WRITE(indcount+100000,'(a19)') "%BLOCK lattice_cart" 
   
           call cell(cellvector, a, b,c,alpha,beta,gamma)
           WRITE(indcount+100000,335) cellvector(1,1), cellvector(2,1), cellvector(3,1)
           WRITE(indcount+100000,335) cellvector(1,2), cellvector(2,2), cellvector(3,2)
           WRITE(indcount+100000,335) cellvector(1,3), cellvector(2,3), cellvector(3,3)
   
           WRITE(indcount+100000,'(a22)') "%ENDBLOCK lattice_cart"
           WRITE(indcount+100000,'(a21)') "%BLOCK positions_frac"
   
   
           do at=1,nat
   		sp=spat(at)
   		if (sp.ne.sptarget) then
                       WRITE(indcount+100000,337) symbol(sp),coords(at,1),coords(at,2),coords(at,3)
   		else
   		    att=at-atini+1
                       call member(nsubs,newconf,att,ifound)
   		    if (ifound==1) then
   		         WRITE (indcount+100000,337) newsymbol(1),coords(at,1),coords(at,2),coords(at,3)
                       else
   		         WRITE (indcount+100000,337) newsymbol(2),coords(at,1),coords(at,2),coords(at,3)
   		    endif
   		endif
 337       format(a3,3(f11.7,2x))
           enddo
   
           WRITE(indcount+100000,'(a24)') "%ENDBLOCK positions_frac"
   
           if (indcount>99999) then
              WRITE (indcount+100000,*) "Error, too many configurations (>99999)! Calculation files not written!"
           endif
   

           do l=1,NLINEAMAX
                READ (62,332,end=499) linea
                WRITE(indcount+100000,333) linea
           enddo
 499    continue
           rewind(62)
   
           if (indcount>99999) then
           WRITE (indcount+100000,*) "Error, too many configurations (>99999)! Calculation files not written!"
           endif

           CLOSE (UNIT=indcount+100000)

        enddo

        CLOSE (UNIT=62)


	END SELECT


!!!!!!!Deallocating arrays
       DEALLOCATE(newconf)
       DEALLOCATE(degen)
       DEALLOCATE(indconf)

       CLOSE (43)

!!!!!!!Reporting the end
       WRITE (*,*) "Done!!!" 
       WRITE (*,*) ""  
       WRITE (*,*) ""  

    END PROGRAM genersod 


