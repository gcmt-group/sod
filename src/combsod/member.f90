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

        SUBROUTINE member(vlength,ivector,iele,f)
        IMPLICIT NONE
        INTEGER,PARAMETER :: NSUBSMAX=20
        INTEGER,DIMENSION(1:NSUBSMAX),INTENT(IN)   :: ivector
        INTEGER,INTENT(IN) :: iele,vlength
        INTEGER,INTENT(OUT) :: f
        INTEGER :: i

! This subrutine checks if iele is a member of ivector

        f=0
        do i=1,vlength
           if  (ivector(i).eq.iele) f=1
        enddo

        RETURN
        END

