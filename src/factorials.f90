FUNCTION pbinomial(nsubs,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, npos 
    INTEGER (kind=8) :: combinations 
    REAL*8 :: x, pbinomial
    pbinomial = combinations(nsubs,npos) * (x**nsubs) * ((1-x)**(npos-nsubs))
END FUNCTION pbinomial


FUNCTION combinations(nsubs, npos)
    IMPLICIT NONE
    INTEGER :: nsubs, npos, k, p
    INTEGER (kind = 8)::  factorial, combinations, ratio

    if (nsubs.le.npos/2) then
         k=nsubs
    else
         k=npos-nsubs
    endif
    ratio=1
    do p=npos-k+1, npos
       ratio=ratio*p
    enddo
    combinations = ratio / factorial(k) 
END FUNCTION combinations


RECURSIVE FUNCTION factorial(n) RESULT (aux)
    IMPLICIT NONE 
    INTEGER, intent(in) :: n
    INTEGER (kind=8):: aux
    if( n == 0 ) then
        aux = 1
    else
        aux = n * factorial(n - 1)
    end if
END FUNCTION

