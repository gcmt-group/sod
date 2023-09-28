FUNCTION momentaA0(nsubsmax,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, nsubsmax, npos
    REAL*8 :: x, pbinomial, momentaA0
    momentaA0 = 0.0
    do nsubs = nsubsmax+1,npos
       momentaA0 = momentaA0 + pbinomial(nsubs,npos,x) 
    enddo
END FUNCTION momentaA0


FUNCTION momentaA1(nsubsmax,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, nsubsmax, npos
    REAL*8 :: x, pbinomial, momentaA1
    momentaA1 = 0.0
    do nsubs = nsubsmax+1,npos
       momentaA1 = momentaA1 + pbinomial(nsubs,npos,x) * nsubs
    enddo
END FUNCTION momentaA1


FUNCTION momentaA2(nsubsmax,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, nsubsmax, npos
    REAL*8 :: x, pbinomial, momentaA2
    momentaA2 = 0.0
    do nsubs = nsubsmax+1,npos
       momentaA2 = momentaA2 + pbinomial(nsubs,npos,x) * (nsubs**2)
    enddo
END FUNCTION momentaA2


FUNCTION momentaB0(nsubsmin,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, nsubsmin, npos
    REAL*8 :: x, pbinomial, momentaB0
    momentaB0 = 0.0
    do nsubs = 0,nsubsmin-1
       momentaB0 = momentaB0 + pbinomial(nsubs,npos,x)
    enddo
END FUNCTION momentaB0


FUNCTION momentaB1(nsubsmin,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, nsubsmin, npos
    REAL*8 :: x, pbinomial, momentaB1
    momentaB1 = 0.0
    do nsubs = 0,nsubsmin-1
       momentaB1 = momentaB1 + pbinomial(nsubs,npos,x) * nsubs
    enddo
END FUNCTION momentaB1


FUNCTION momentaB2(nsubsmin,npos,x)
    IMPLICIT NONE
    INTEGER :: nsubs, nsubsmin, npos
    REAL*8 :: x, pbinomial, momentaB2
    momentaB2 = 0.0
    do nsubs = 0,nsubsmin-1
       momentaB2 = momentaB2 + pbinomial(nsubs,npos,x) * (nsubs**2)
    enddo
END FUNCTION momentaB2


