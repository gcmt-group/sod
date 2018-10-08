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

        SUBROUTINE Bubble(A, N)

	IMPLICIT NONE
	INTEGER,PARAMETER :: NATMAX=1000
	INTEGER :: I,J,ITemp,N
	INTEGER,DIMENSION(1:NATMAX)   :: A

!
!       Make N-1 passes through the array.
!       On pass i, "bubble" the next smallest element
!       up from the end of the array to position i.
!
        DO I = 1, N-1
        DO J = N, I+1, -1
            IF (A(J) < A(J-1)) THEN
              CALL Swap(A(J), A(J-1))
            ENDIF
          END DO 
        END DO 
        RETURN
        END

        SUBROUTINE Swap(I, J)
!       Exchanges the integers I and J.
!
        ITemp = I
        I = J
        J = ITemp
        RETURN
        END

