
!    Copyright (c) 2022 Ricardo Grau-Crespo, Said Hamad
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

    PROGRAM combsod 
       IMPLICIT NONE

       INTEGER,PARAMETER :: NSPMAX=10 , NATMAX=10000, NOPMAX=10000, NCELLMAX=1000
       INTEGER,PARAMETER :: NLINEAMAX=200 
       REAL,PARAMETER :: tol0=0.001, kB=8.61734E-5

       INTEGER :: i,j,k,l,xx,yy,suma,t,ina,inb,inc,elei, sub
       INTEGER :: op1, nop1, op, nop, opsc, nopsc,ifound,op1new,nop1new,aux,opc,nopc
       INTEGER :: sp, nsp, spmap, nspmap, cumnatsp, cumnatspmap, sptarget, sptargetmap,ssp
       INTEGER :: at0, nat0, at1, nat1, nat1r, at, nat , atmap, natmap, at1r, at1i,attmp,att,mapno,FILER, MAPPER
       INTEGER :: na,nb,nc,nsubs,nsubsmap,atini,atfin,atinimap,atfinmap
       INTEGER :: pos,npos
       INTEGER (kind=8):: ntc, count, nic, equivcount, iequiv, indcount, combinations
       LOGICAL :: found,foundnoind,mores
       INTEGER,DIMENSION(:), ALLOCATABLE:: newconf
       INTEGER,DIMENSION(2)   :: newshell, newshellmap
       INTEGER,DIMENSION(NOPMAX)   :: op1good
       REAL,DIMENSION(NOPMAX,3,3) :: mgroup1, mgroup, mgroup1new
       REAL,DIMENSION(NOPMAX,3)   :: vgroup1, vgroup,vgroup1sca,vgroup1scb,vgroup1scc,vgroup1new
       INTEGER,DIMENSION(NOPMAX,NATMAX)   :: fulleqmatrix,eqmatrixtarget
       INTEGER,DIMENSION(NATMAX) :: spat0, spat1, spat, spatmap, spat1r
       REAL,DIMENSION(NATMAX,3) :: coords0, coords1, coords, coords1r, result, coordsmap
       INTEGER,DIMENSION(NATMAX)   :: as,test
       INTEGER,DIMENSION(:), ALLOCATABLE:: degen
       INTEGER,DIMENSION(NSPMAX) :: natsp0, natsp1, natsp, natspmap, snatsp, ishell, ishellmap
       REAL,DIMENSION(3) :: coordstemp, vcorr
       REAL,DIMENSION(3,3) :: cellvector
       REAL,DIMENSION(NCELLMAX,3)   :: vt
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: conf,indconf,indconfmap,equivconf
       INTEGER,DIMENSION(:,:), ALLOCATABLE:: mapping
       REAL :: a1,b1,c1,alpha,beta,gamma,a,b,c,cc, amap, bmap, cmap, alphamap, betamap, gammamap
       REAL :: prod, x, maxentropy, ientropy, perc
       CHARACTER,DIMENSION(NSPMAX) :: symbol*3, symbolmap*3, ssymbol*3
       CHARACTER,DIMENSION(2) :: newsymbol*3, newsymbolmap*3
       CHARACTER :: linea*85, title*15, tmptitle*15, runtitle*40
       CHARACTER :: sindcount*5, snumber*4, tmpstr*50 


! Input files

       OPEN (UNIT= 9,FILE="INSOD")
       OPEN (UNIT=12,FILE="SGO")

! Output files

       OPEN (UNIT=25,FILE="coordinates.xyz")
       OPEN (UNIT=26,FILE="EQMATRIX")
       OPEN (UNIT=30,FILE="OUTSOD")
       OPEN (UNIT=31,FILE="SUPERCELL")
       OPEN (UNIT=43,FILE="filer")
       OPEN (UNIT=46,FILE="OPERATORS")
       OPEN (UNIT=47,FILE="cSGO")

!
! DEFINITION OF VARIABLES:
!
! op                  Index for the operators in the supercell
! op1                 Index for the operators in the unit cell
! nop                 Total number of operators in the supercell
! nop1                Total number of operators in the unit cell
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
       WRITE (*,*) "         SOD (Site Occupancy Disorder) version 0.52  " 
       WRITE (*,*) " " 
       WRITE (*,*) "         Authors: R. Grau-Crespo and S. Hamad                                   " 
       WRITE (*,*) " " 
       WRITE (*,*) "         Contact:  <r.grau-crespo@reading.ac.uk> " 
       WRITE (*,*) "**************************************************************************** " 
       WRITE (*,*) " " 
       WRITE (*,*) " " 
       WRITE (*,*) " " 
       WRITE (*,*) "Reading input files..." 
       WRITE (*,*) " " 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Reading the input file with the uc space group information: SGO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       READ (12,*)
       READ (12,*) op1
       do while(op1>0)
        do i=1,3
          READ (12,*)(mgroup1(op1,i,j),j=1,3),vgroup1(op1,i)
         enddo
        nop1=op1
        READ (12,*) op1
       enddo


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


!cccccccccccccccccccccccccccccccccccc
! Generating spat0 array
!cccccccccccccccccccccccccccccccccccc

       do at0=1, natsp0(1)
          spat0(at0)=1
          cumnatsp=natsp0(1)
       enddo 
       do sp=2, nsp
          do at0=cumnatsp+1, cumnatsp+natsp0(sp)
             spat0(at0)=sp
          enddo
          cumnatsp=cumnatsp+natsp0(sp) 
       enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      First create the redundant unit cell from the asymmetric unit 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
       WRITE (*,*) "" 
       WRITE (*,*) "Generating the supercell..." 
       WRITE (*,*) "" 

       at1r=0
       do at0=1,nat0
          do op1=1,nop1
             at1r=at1r+1
             coords1r(at1r,1:3)= MATMUL(mgroup1(op1,1:3,1:3),coords0(at0,1:3))+vgroup1(op1,1:3)
             coords1r(at1r,1)=cc(coords1r(at1r,1))
             coords1r(at1r,2)=cc(coords1r(at1r,2))
             coords1r(at1r,3)=cc(coords1r(at1r,3))
             spat1r(at1r)=spat0(at0)
          enddo
       enddo
       nat1r=at1r

 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Get rid of the redundant atoms, to create the unit cell 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	coords1(1,:)=coords1r(1,:)
	at1r=1
	at1=1
	coords1(1,:)=coords1r(1,:)
	spat1(1)=spat1r(1)
	do at1r=2,nat1r
           found=.FALSE.
	   do at1i=1,at1
              prod=DOT_PRODUCT(coords1r(at1r,:)-coords1(at1i,:),&
                            coords1r(at1r,:)-coords1(at1i,:)) 
	      if (prod.le.tol0) found=.TRUE.
	   enddo
           if (found.eqv..FALSE.) then    
	       at1=at1+1
	       coords1(at1,:)=coords1r(at1r,:)
		spat1(at1)=spat1r(at1r)
	   endif
        enddo  
        nat1=at1

	do sp=1,nsp
	   natsp1(sp)=0
	enddo
	do at1=1,nat1
	   natsp1(spat1(at1))=natsp1(spat1(at1))+1
	enddo

	natsp(:)=na*nb*nc*natsp1(:)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create the traslation vectors of the supercell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	nopsc=nop1
	t=0
	do ina=0,na-1
	  do inb=0,nb-1
	    do inc=0,nc-1
		t=t+1
		vt(t,1)=real(ina)/real(na)
		vt(t,2)=real(inb)/real(nb)
		vt(t,3)=real(inc)/real(nc)
	    enddo 
	  enddo
	enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create the traslation vectors of the supercell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do op1=1,nop1
	   op1good(op1)=1
	enddo

	if (na.eq.nb) then 
	   if (na.eq.nc) then 
	      goto 255
	      else
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,1,3).ne.0).OR.(mgroup1(op1,3,1).ne.0).OR.&
		        (mgroup1(op1,2,3).ne.0).OR.(mgroup1(op1,3,2).ne.0)) then
	            op1good(op1)=0
		    endif
                 enddo
              enddo
	   endif
	else
	   if (na.eq.nc) then 
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,1,2).ne.0).OR.(mgroup1(op1,2,1).ne.0).OR.&
		        (mgroup1(op1,3,2).ne.0).OR.(mgroup1(op1,2,3).ne.0)) then
	            op1good(op1)=0
	            endif
                 enddo
              enddo
	   endif
	   if (nb.eq.nc) then 
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,2,1).ne.0).OR.(mgroup1(op1,1,2).ne.0).OR.&
		        (mgroup1(op1,3,1).ne.0).OR.(mgroup1(op1,1,3).ne.0)) then
	            op1good(op1)=0
	            endif
                 enddo
              enddo
	   endif
	   if ((nb.ne.nc).AND.(na.ne.nc)) then 
              do op1=1,nop1
                 do i=1,3
                    if ((mgroup1(op1,2,1).ne.0).OR.(mgroup1(op1,1,2).ne.0).OR.&
		        (mgroup1(op1,3,1).ne.0).OR.(mgroup1(op1,1,3).ne.0).OR.&
		        (mgroup1(op1,3,2).ne.0).OR.(mgroup1(op1,2,3).ne.0)) then
	            op1good(op1)=0
	            endif
                 enddo
              enddo
	   endif
	endif      
 255       continue

	op1new=0
        do op1=1,nop1
	   if (op1good(op1).eq.1) then
 	       op1new=op1new+1
	       mgroup1new(op1new,:,:)=mgroup1(op1,:,:)
	       vgroup1new(op1new,:)=vgroup1(op1,:)
	   endif
	enddo
	nop1new=op1new

	nop1=nop1new
	do op1=1,nop1
           mgroup1(op1,:,:)=mgroup1new(op1,:,:)
           vgroup1(op1,:)=vgroup1new(op1,:)
	enddo



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Generate the supercell coordinates, by applying these vectors
!      to the unit cell. Also calculate the supercell parameters.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	a=na*a1
	b=nb*b1
	c=nc*c1


	at=0
	do at1=1,nat1
	   do t=1,na*nb*nc
		at=at+1
		coords(at,1)=vt(t,1)+coords1(at1,1)/na
		coords(at,2)=vt(t,2)+coords1(at1,2)/nb
		coords(at,3)=vt(t,3)+coords1(at1,3)/nc
                spat(at)=spat1(at1)
	   enddo
	enddo
	nat=at

	at=0
	do at1=1,nat1
	   do t=1,na*nb*nc
		at=at+1
 100    format(3(f11.7,2x))
	   enddo
	enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Write a file with the fractional coordinates of the supercell
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        WRITE(31,111) a,b,c,alpha,beta,gamma 
 111    format(6(f10.4,2x))
        WRITE(31,*) natsp(1:nsp)
        do at=1,nat
                 WRITE(31,110) symbol(spat(at)),coords(at,1),coords(at,2),coords(at,3)
 110    format(a3,2x,3(f11.7,2x))
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate the supercell operators
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	op=0
	do op1=1,nop1
	   do t=1,na*nb*nc
		op=op+1
		mgroup(op,:,:)=mgroup1(op1,:,:)		
		vgroup(op,1)=vgroup1(op1,1)/na+vt(t,1)
		vgroup(op,2)=vgroup1(op1,2)/nb+vt(t,2)
		vgroup(op,3)=vgroup1(op1,3)/nc+vt(t,3)
!		WRITE(*,*) op,op1,t
	   enddo
	enddo

	nop=op

	WRITE(46,*) nop
	do op=1,nop
	   WRITE(46,*) "Operator number ",op
	   do i=1,3
	   WRITE(46,*) mgroup(op,i,1:3),vgroup(op,i)
 163   format(4(f10.8))
	   enddo
	enddo
       
       WRITE (*,*) "       Composition of the original supercell (parent structure):" 
       do sp=1,nsp
          WRITE (*,*) "                                                         ", symbol(sp),natsp(sp)
       enddo     
       WRITE (*,*) " "
       WRITE (*,*) "       Number of symmetry operators in the supercell:       ", nop 
       WRITE (*,*) " " 
       WRITE (*,*) "       Composition of the substituted supercell:" 
       do sp=1,nsp
          if (sp==sptarget) then
             WRITE (*,*) "                                                         ", newsymbol(1),nsubs
             WRITE (*,*) "                                                         ", newsymbol(2),natsp(sp)-nsubs
          else   
             WRITE (*,*) "                                                         ", symbol(sp),natsp(sp)
          endif
       enddo

       WRITE (*,*) ""
       WRITE (*,*) ""
       WRITE (*,*) "Generating the complete configurational space..." 
       WRITE (*,*) " " 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Create Full Equivalence Matrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	coordstemp(:)=MATMUL(mgroup(1,:,:),coords(1,:))+vgroup(1,:)

	do op=1,nop
	   do at=1,nat
	        coordstemp(:)=MATMUL(mgroup(op,:,:),coords(at,:))+vgroup(op,:)
                coordstemp(1)=cc(coordstemp(1))
                coordstemp(2)=cc(coordstemp(2))
                coordstemp(3)=cc(coordstemp(3))
	        found=.FALSE.	
		attmp=0
		do while((found.eqv..FALSE.).AND.(attmp.lt.nat))
		    attmp=attmp+1
                    if((abs(coordstemp(1)-coords(attmp,1)).lt.tol0).AND.&
                       (abs(coordstemp(2)-coords(attmp,2)).lt.tol0).AND.&
                       (abs(coordstemp(3)-coords(attmp,3)).lt.tol0))& 
                       found=.TRUE.
	    enddo
		If( (attmp.eq.nat).And.(.Not.found) ) Then 
             Write(*,*) "Error!!! Operator",& 
             	        op,"applied on atom",at,"does not produce another atom in the list!!!"
             attmp=0
        Endif

        fulleqmatrix(op,at)=attmp

	   enddo
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Obtain the Equivalence Matrix for the target species from the Full Equivalence Matrix
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	npos=natsp(sptarget)
        if (npos==nsubs) then
            write(*,*) "Ilegal number of substitutions"
            STOP
        endif
	do op=1,nop
	att=0
	do at=atini,atfin
	   att=att+1
           eqmatrixtarget(op,att)=fulleqmatrix(op,at)-atini+1
	enddo
	enddo


  WRITE(26,*) nop,npos
  do op=1, nop
      WRITE(26,*) (eqmatrixtarget(op,pos), pos=1,npos)
  enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Generate the list of configurations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
        ntc=combinations(nsubs,npos)
        maxentropy= kB*LOG(REAL(ntc))	
        WRITE(*,*) " " 
        WRITE(*,*) "       Total number of configurations in the supercell:     ", ntc 
        WRITE(*,*) " " 
        WRITE(*,*) "       Maximum entropy for this composition and supercell:", maxentropy, " eV/K"
        WRITE(*,*) " " 


        x=REAL(nsubs)/REAL(npos)
        WRITE(*,*) "       Fraction of substituted sites:                 x = ", x
        WRITE(*,*) " " 
        ientropy=-npos*kB*(x*log(x)+(1-x)*log(1-x))
        WRITE(*,*) "       Ideal entropy (per cell) for this composition:     ", ientropy, " eV/K"
        WRITE(*,*) " " 
!        ntcmax=exp(ientropy/kB)
        
!!!!!!!!Allocating array sizes

        ALLOCATE(degen(1:ntc))
        ALLOCATE(newconf(1:nsubs))
        ALLOCATE(conf(1:ntc,1:nsubs))
        ALLOCATE(indconf(1:ntc,1:nsubs))
        ALLOCATE(equivconf(1:ntc,1:nsubs))

        mores=.false.
        as(:)=1
        count=1
        call ksubset(npos,nsubs,as,mores)
        conf(1,1:nsubs)=as(1:nsubs)
        do while(mores)
           call ksubset(npos,nsubs,as,mores)
           count=count+1
           conf(count,1:nsubs)=as(1:nsubs)
        enddo

        if (count.ne.ntc) then
                write(*,*) "Error in KSUBSET subroutine"
                STOP
        endif


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Generate the list of independent configurations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        WRITE (*,*) " " 
        WRITE (*,*) "Finding the inequivalent configurations..." 
        WRITE (*,*) " " 
        WRITE (*,*) "       Found    Completion " 
        WRITE (*,*) "       =====    ========== " 

! With this loop we get a function that returns the configuration equivalent, by the operator op,
! to the configuration count.

! For the first configuration (count=1)	
        indconf(:,:)=0
        equivconf(:,:)=0
        newconf(:)=0
	degen(:)=1

        count=1
        equivcount=1
        indcount=1
        WRITE(47,*) "List of operators for configuration: ",indcount

        indconf(indcount,1:nsubs)=conf(count,1:nsubs)
	opc=1
        WRITE (47,*) opc
        do i=1,3
          WRITE (47,*)(mgroup(1,i,j),j=1,3),vgroup(1,i)
        enddo
        equivconf(equivcount,1:nsubs)=conf(count,1:nsubs)
        do op=2,nop
                do i=1,nsubs
                    elei=conf(count,i)                 ! elei is the element i of the configuration number count
                    newconf(i)=eqmatrixtarget(op,elei) ! newconf is the configuration obteined by applying the 
					               ! operator op to the configuration number conf(count,i)
                enddo
                call bubble(newconf,nsubs)
                found=.FALSE.
                prod=DOT_PRODUCT(newconf(1:nsubs)-conf(count,1:nsubs),&
                                 newconf(1:nsubs)-conf(count,1:nsubs))
                if (prod.le.tol0) then
	            found=.TRUE.
		    opc=opc+1           
                    WRITE (47,*) opc
         	    do i=1,3
                       WRITE (47,*)(mgroup(op,i,j),j=1,3),vgroup(op,i)
         	    enddo
                else 
                    do iequiv=1,equivcount
                       prod=DOT_PRODUCT(newconf(1:nsubs)-equivconf(iequiv,1:nsubs),&
                                        newconf(1:nsubs)-equivconf(iequiv,1:nsubs))
                       if (prod.le.tol0) found=.TRUE.
                    enddo
	        endif	
                if (found.eqv..FALSE.) then
                    equivcount=equivcount+1
                    equivconf(equivcount,1:nsubs)=newconf(1:nsubs)
                    degen(indcount)=degen(indcount)+1
                endif
        enddo
        WRITE (47,*) 0

! For the rest of configurations (count>1)	
        perc=0.0
        do while(equivcount.lt.ntc)
           count=count+1
           foundnoind=.FALSE.
           do iequiv=1,equivcount
                 prod=DOT_PRODUCT(conf(count,1:nsubs)-equivconf(iequiv,1:nsubs),&
                                  conf(count,1:nsubs)-equivconf(iequiv,1:nsubs))
                 if (prod.le.tol0) foundnoind=.TRUE.
           enddo
           if (foundnoind.eqv..FALSE.) then
               indcount=indcount+1
               WRITE(47,*) "List of operators for configuration: ",indcount
               indconf(indcount,1:nsubs)=conf(count,1:nsubs)
               !!!!!!!!!!!    UPDATING EQUIVCONF  !!!!!!!!!!!!!!!
               equivcount=equivcount+1
               equivconf(equivcount,1:nsubs)=conf(count,1:nsubs)
	       opc=1
               WRITE (47,*) opc
               do i=1,3
                 WRITE (47,*)(mgroup(1,i,j),j=1,3),vgroup(1,i)
               enddo
               do op=2,nop
                   do i=1,nsubs
                       elei=conf(count,i)
                       newconf(i)=eqmatrixtarget(op,elei)
                   enddo
                   call bubble(newconf,nsubs)
                   found=.FALSE.
                   prod=DOT_PRODUCT(newconf(1:nsubs)-conf(count,1:nsubs),&
                                    newconf(1:nsubs)-conf(count,1:nsubs))
                   if (prod.le.tol0) then
                       found=.TRUE.
		       opc=opc+1           
                       WRITE (47,*) opc
         	       do i=1,3
                          WRITE (47,*)(mgroup(op,i,j),j=1,3),vgroup(op,i)
         	       enddo
                   else
	               do iequiv=1,equivcount
                           prod=DOT_PRODUCT(newconf(1:nsubs)-equivconf(iequiv,1:nsubs),&
                                            newconf(1:nsubs)-equivconf(iequiv,1:nsubs))
                           if (prod.le.tol0) then
			     found=.TRUE.
              		     goto 200
                           endif
                       enddo
	           endif
                   if (found.eqv..FALSE.) then
                       equivcount=equivcount+1
                       equivconf(equivcount,1:nsubs)=newconf(1:nsubs)
                       degen(indcount)=degen(indcount)+1
                   endif
 200               continue
               enddo
               WRITE (47,*) 0
               if ((100.0*equivcount/ntc - perc > 5.0).OR.(equivcount==ntc)) then
                  perc=100.0*equivcount/ntc
                  WRITE(*,'(4x,i6,7x,f5.1,a2)') indcount, perc, "% "
               endif 
           endif
	enddo

	nic=indcount

			
	WRITE(*,*) " "
	WRITE(*,*) "       Number of inequivalent configurations:               ",nic
	WRITE(*,*) " "

  WRITE(30,*) nsubs, " substitutions in ", npos, "sites" 
  WRITE(30,*) nic, " configurations" 
	do indcount=1,nic
	WRITE(30,330) indcount,degen(indcount),indconf(indcount,1:nsubs)
 330    format(i6,1x,i6,30(1x,i4))
	enddo

!!!!!!!!End of search for independent configurations



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!This part of the program is only accessed if mapping is required
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       if (MAPPER>0) then
 
           ALLOCATE(mapping(1:npos,1:MAPPER))
           ALLOCATE(indconfmap(1:ntc,1:nsubs*MAPPER))


!!!!!!!!!!!Opening and reading MAPTO file
           
           OPEN (UNIT=13,FILE="MAPTO")
                                 
           READ (13,*)
           READ (13,*) 
           READ (13,*)
           READ (13,*)
           READ (13,*) amap,bmap,cmap,alphamap,betamap,gammamap
           READ (13,*)
           READ (13,*)
           READ (13,*) nspmap
           if (nspmap.ne.nsp) WRITE(*,*) "Warning! Different number of species in original and new structures."
           READ (13,*)
           READ (13,*)
           READ (13,*)(symbolmap(spmap),spmap=1,nspmap)
           READ (13,*)
           READ (13,*)
           READ (13,*)(natspmap(spmap),spmap=1,nspmap)
           READ (13,*)
           READ (13,*)

           natmap=0
           do spmap=1,nspmap
              natmap = natmap + natspmap(spmap)
           enddo

           do atmap=1,natmap
             READ (13,*)(coordsmap(atmap,i),i=1,3)
           enddo


           READ (13,*)
           READ (13,*)
           READ (13,*) sptargetmap
           READ (13,*)
           READ (13,*)
           READ (13,*)
           READ (13,*) (newsymbolmap(i),i=1,2)
           READ (13,*)
           READ (13,*)
           READ (13,*)
           READ (13,*) (ishellmap(sp),sp=1,nsp)
           READ (13,*)
           READ (13,*) (newshellmap(i),i=1,2)
           READ (13,*)
           READ (13,*)
           READ (13,*)
          
           do pos=1,npos
             READ (13,*)(mapping(pos,mapno), mapno=1,MAPPER)
           enddo


!!!!!!!!obtaining some descriptors of the new structure

           !spatmap array
           do atmap=1, natspmap(1)
              spatmap(atmap)=1
           enddo
           cumnatspmap=natspmap(1)
           do spmap=2, nspmap
              do atmap=cumnatspmap+1, cumnatspmap+natspmap(spmap)
                 spatmap(atmap)=spmap
              enddo
              cumnatspmap=cumnatspmap+natspmap(spmap)
           enddo


           !atinimap and atfinmap indexes
           if (sptargetmap.eq.1) then
              atinimap=1
           else
              atinimap=1
              do spmap=1,sptargetmap-1
                 atinimap=atinimap+natspmap(spmap)
              enddo
           endif
           atfinmap=atinimap+natspmap(sptargetmap)-1


!!!!!!!!!!Obtaining the new list of independent configurations

           do indcount=1,nic
              do mapno=1,MAPPER
                 do sub=1,nsubs
                    indconfmap(indcount,(mapno-1)*nsubs+sub)=mapping(indconf(indcount,sub),mapno)
                 enddo
              enddo
           enddo


!!!!!!!!!!Giving the new values to the old variables
           
           a=amap
           b=bmap
           c=cmap
           alpha=alphamap
           beta=betamap
           gamma=gammamap
           nat=natmap
           nsp=nspmap
           natsp(1:nsp)=natspmap(1:nsp)
           spat(1:nat)=spatmap(1:nat) 
           sptarget=sptargetmap
           coords(1:nat,1:3)=coordsmap(1:nat,1:3)           
           symbol(1:nsp)=symbolmap(1:nsp)
           newsymbol(1:2)=newsymbolmap(1:2)
           ishell(1:nsp)=ishellmap(1:nsp)
           newshell(1:2)=newshellmap(1:2)
           atini=atinimap
           atfin=atfinmap
           
 
!!!!!!!!!!!and only now, not before          
           nsubs=MAPPER*nsubs
           npos=MAPPER*npos


!!!!!!!!!!!Reallocating arrays needed for printing files
           DEALLOCATE(newconf)
           DEALLOCATE(indconf)
           ALLOCATE(newconf(1:nsubs))
           ALLOCATE(indconf(1:ntc,1:nsubs))

           indconf(:,:)=indconfmap(:,:) 

!!!!!!!!!!!Deallocating arrays no longer needed
           DEALLOCATE(mapping)
           DEALLOCATE(indconfmap)


       endif


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

        if (indcount<10) then
        WRITE(UNIT=snumber, FMT='(I1)') indcount
        WRITE (indcount+100000,334) "output xtl c0000",snumber
        WRITE (indcount+100000,334) "output arc c0000",snumber
        endif
        if ((indcount>9).AND.(indcount<100)) then
        WRITE(UNIT=snumber, FMT='(I2)') indcount
        WRITE (indcount+100000,334) "output xtl  c000",snumber
        WRITE (indcount+100000,334) "output arc  c000",snumber
        endif
        if ((indcount>99).AND.(indcount<1000)) then
        WRITE(UNIT=snumber, FMT='(I3)') indcount
        WRITE (indcount+100000,334) "output xtl   c00",snumber
        WRITE (indcount+100000,334) "output arc   c00",snumber
        endif
        if ((indcount>999).AND.(indcount<10000)) then
        WRITE(UNIT=snumber, FMT='(I4)') indcount
        WRITE (indcount+100000,334) "output xtl    c0",snumber
        WRITE (indcount+100000,334) "output arc    c0",snumber
        endif
        if ((indcount>9999).AND.(indcount<100000)) then
        WRITE(UNIT=snumber, FMT='(I4)') indcount
        WRITE (indcount+100000,334) "output xtl    c",snumber
        WRITE (indcount+100000,334) "output arc    c",snumber
        endif

        if (indcount>99999) then
	WRITE (indcount+100000,*) "Error, too many configurations (>99999)! Calculation files not written!"
	endif
 334    format(a15,a4)

        CLOSE (UNIT=indcount+100000)

        enddo

        CLOSE (UNIT=41)
        CLOSE (UNIT=42)

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
        CLOSE (UNIT=51)
        CLOSE (UNIT=52)

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



!!!!!!!Reporting the end
       WRITE (*,*) "Done!!!" 
       WRITE (*,*) ""  
       WRITE (*,*) ""  

    END PROGRAM combsod 


