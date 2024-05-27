      PROGRAM AMCABC

      !Variable declaration
      IMPLICIT real*8 (A-H,P-V,X-Z)
      IMPLICIT integer (I-N)
      IMPLICIT logical (O)
      IMPLICIT character*150 (W)

      PARAMETER (npA=3, npB=0, npC=0)
      PARAMETER (mA=4, mB=1, mC=1)
      PARAMETER (iDir=1, iScr=66)
      PARAMETER (wFolder = 'VMC_s4_e12_ags1')

      PARAMETER (hbar = 1.054571817D-34)
      PARAMETER (bk = 1.380649D-23)
      PARAMETER (u = 1.660539066D-27)

      PARAMETER (n = 1)
      DIMENSION Ck(3), wTwf(6), sig(6), eps4(6)
      integer,DIMENSION (:,:), ALLOCATABLE :: ijp
      include 'include/TWF_dim.for'

      real*8,DIMENSION (:), ALLOCATABLE :: SEk, SEp, Ep, ePnew, ePold
      real*8,DIMENSION (:,:), ALLOCATABLE :: V, Ek
      real*8,DIMENSION (:,:), ALLOCATABLE :: Aux, r1, r2
      real*8,DIMENSION (:,:,:), ALLOCATABLE :: F, rn, rp, dr

      !include 'include/VMC_P/R_1D_0.4.for' !Declare

************************************************************************
* 0 **    Setting initial parameters                              ******
c OVDJE JE KONVERZIJA JEDINICA
      !Ck = -h_bar**2/(2m)  [mK A**2]
      Ck(1)=-1.d0*hbar**2/(2.d0*mA*u)/(1.d-20*bk)*1.d3
c 1.d3 znači 1 puza deset na treću
      Ck(2)=-1.d0*hbar**2/(2.d0*mB*u)/(1.d-20*bk)*1.d3
      Ck(3)=-1.d0*hbar**2/(2.d0*mC*u)/(1.d-20*bk)*1.d3

      npAB = npA + npB
      np  = npAB + npC
      include 'include/DO_dir_set.for'! npA(.npB(.npC))
      include 'include/VMC_set.for'
      include 'include/E_Ep_set.for'
      include 'include/TWF.for'

      DO iRun = 1, nRun !Restart enabled
      include 'include/VMC_start.for'

************************************************************************
* 1 **  Starting positons and energies                            ******
      do iw = 1, nw
        do ip = 1, np
          do k=1,2
            rn(k,ip,n)=rp(k,ip,iw)
          enddo
        enddo
        include 'include/DO_calc_R1_R2.for'
        include 'include/TWF_sum.for'
        ePold(iw) = ePnew(iw)
        !Calculating energy
        include 'include/E_Ep.for'
        include 'include/E_Ek.for'
      enddo
      !Reseting variables to initial values
      SbEk = 0.d0
      SbEp = 0.d0
      SbE2 = 0.d0
      nAcc = 0.d0
      SbAcc = 0.d0
      Ecut = 3.0d2
      !include 'include/VMC_P/R_1D_1.4.for' !Initial distances
      do ib = 1, nb
       SswEp = 0.d0
       SswEk = 0.d0
       ob = ib .GT. nbSkip
       ibEff = ib - nbSkip
       do is = 1, ns
        do iw = 1, nw

************************************************************************
* 2 **    Comparing probabilities of old and new position         ******
          !Moving particle (index n is auxiliary index needed in DMC)
          do ip = 1, np
            do k=1,2
              rn(k,ip,n)=rp(k,ip,iw)+step*(ran1(nrg)-0.5d0)*2.d0
            enddo
          enddo
          include 'include/DO_calc_R1_R2.for'
          include 'include/TWF_sum.for'
          tmp = ePnew(iw) - ePold(iw)
          if(tmp .GE. 0.d0) then
            oAcc = .TRUE.
          elseif( ran1(nrg) .LE. DEXP(2.d0*tmp)) then
            oAcc = .TRUE.
          else
            oAcc = .FALSE.
          endif

************************************************************************
* 3 **    Calculation of local energy and P for new position     ******
          if(oAcc) then
            nAcc = nAcc + 1       !Calculating acceptance in block
            include 'include/E_Ep.for'
            include 'include/E_Ek.for'
            !New position becomes old one
            do ip = 1, np
              do k = 1, 2
                rp(k,ip,iw)=rn(k,ip,n)
              enddo
            enddo
            ePold(iw) = ePnew(iw)
          !include 'include/VMC_P/R_1D_2.4.for' !New distances
          endif

************************************************************************
* 4 **    Sumation of loacl energy                                ******
          SswEk = SswEk + SEk(iw)
          SswEp = SswEp + SEp(iw)
          !include 'include/VMC_P/R_1D_3.4.for' !Sorting distributions
        enddo !walkers
       enddo !steps in block
       !Changig step lenght according acceptance in a block
       Acc = dfloat(nAcc)/dfloat(nw*ns)
       if( Acc .GT. 0.5d0) then
         step = step * 1.05d0
       else
         step = step * 0.95d0
       endif
       if(ob) then
        SbAcc = SbAcc + Acc
        !Averaging energies
        AbEp=SswEp/(nw*ns)
        AbEk=SswEk/(nw*ns)
        SbEp=SbEp+AbEp
        SbEk=SbEk+AbEk
        SbE2=SbE2+(AbEk+AbEp)**2
        include 'include/E_save_VMC.for'


       endif
       nAcc = 0
      enddo !blocks


      include 'include/VMC_end_R.for'
      include 'include/VMC_end_E.for'
      !include 'include/VMC_P/R_1D_4.4.for' !Saving
      !Read saved positions and run again (loadr 1->2)
      loadr = loadr +1
      step = sTemp
      nbSkip = 0
      nb = nbTmp
      ENDDO


      CLOSE(iScr)
      STOP
      END


      !NR random number generator
      include 'include/ran2.for'

