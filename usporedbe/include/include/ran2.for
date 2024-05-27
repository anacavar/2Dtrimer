************************************************************************
**                                                                    **
**  (C) Copr. 1986-92 Numerical Recipes Software *5jw#.               **
** Random number generator                                            **
** which generates nubers from <0,1> distribuded uniformly            **
**                                                                    **
************************************************************************
      real*8 FUNCTION ran1(nrg)
      INTEGER nrg,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,
     &           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     &           IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-15,
     &           RNMX=1.d0-EPS)
      INTEGER nrg2,j,k,iv(NTAB),iy
      SAVE iv,iy,nrg2
      DATA nrg2/123456789/, iv/NTAB*0/, iy/0/
      if (nrg.le.0) then
        nrg=max(-nrg,1)
        nrg2=nrg
        do 11 j=NTAB+8,1,-1
          k=nrg/IQ1
          nrg=IA1*(nrg-k*IQ1)-k*IR1
          if (nrg.lt.0) nrg=nrg+IM1
          if (j.le.NTAB) iv(j)=nrg
 11     continue
        iy=iv(1)
      endif
      k=nrg/IQ1
      nrg=IA1*(nrg-k*IQ1)-k*IR1
      if (nrg.lt.0) nrg=nrg+IM1
      k=nrg2/IQ2
      nrg2=IA2*(nrg2-k*IQ2)-k*IR2
      if (nrg2.lt.0) nrg2=nrg2+IM2
      j=1+iy/NDIV
      iy=iv(j)-nrg2
      iv(j)=nrg
      if(iy.lt.1)iy=iy+IMM1
      ran1=min(AM*iy,RNMX)
      return
      END
