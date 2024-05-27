      PROGRAM DMCABC

      !Variable declaration
      IMPLICIT real*8 (A-H,P-V,X-Z)
      IMPLICIT integer (I-N)
      IMPLICIT logical (O)
      IMPLICIT character*150 (W)

      PARAMETER (npA=3, npB=0, npC=0)
      PARAMETER (mA=4, mB=4, mC=12)
      PARAMETER (iDir=1, iScr=66)
      PARAMETER(iTau=10, powTau=1.d-7)
      PARAMETER(wFolder = 'DMC_s4_e6.0_E_8kw')

      PARAMETER (hbar = 1.054571628D-34)
      PARAMETER (bk = 1.3806504D-23)
      PARAMETER (u = 1.660538782D-27)

      DIMENSION Ck(3), rCk(3), var(3), wTwf(6), sig(6), eps4(6)
      integer,DIMENSION (:), ALLOCATABLE :: nsx
      integer,DIMENSION (:,:), ALLOCATABLE :: ijp
      character*50,DIMENSION (:), ALLOCATABLE :: wSuf
      include 'include/TWF_dim.for'
      real*8,DIMENSION (:), ALLOCATABLE :: SEk, SEp, Ep
      real*8,DIMENSION (:,:), ALLOCATABLE :: V, Ek, Elocal
      real*8,DIMENSION (:,:), ALLOCATABLE :: Aux, r1, r2
      real*8,DIMENSION (:,:,:), ALLOCATABLE :: F, rn, dr
      real*8,DIMENSION (:,:,:,:), ALLOCATABLE :: rp

      !include 'include/pd_pur_R2/pd_pur.0.5.for' !dimension

************************************************************************
* 0 **    Setting initial parameters                              ******
      !Ck = -h_bar**2/(2m)  [mK A**2]
      Ck(1)=-1.d0*hbar**2/(2.d0*mA*u)/(1.d-20*bk)*1.d3
      Ck(2)=-1.d0*hbar**2/(2.d0*mB*u)/(1.d-20*bk)*1.d3
      Ck(3)=-1.d0*hbar**2/(2.d0*mC*u)/(1.d-20*bk)*1.d3

      rCk(1)=1.d0/Ck(1)
      rCk(2)=1.d0/Ck(2)
      rCk(3)=1.d0/Ck(3)

      tau  = iTau*powTau
      var(1) = dsqrt(2.d0*Ck(1)*(-1.d0)*tau)
      var(2) = dsqrt(2.d0*Ck(2)*(-1.d0)*tau)
      var(3) = dsqrt(2.d0*Ck(3)*(-1.d0)*tau)

      npAB = npA + npB
      np  = npAB + npC
      n=1
      include 'include/DO_dir_set.for'
      include 'include/DMC_set.for'
      !include 'include/pd_pur_R2/pd_pur.1.5.for' !allocating 
      include 'include/E_Ep_set.for'
      include 'include/TWF.for'

************************************************************************
* 1 **    For each walker read its position and loacl energy      *****
      include 'include/DMC_start.for'
      include 'include/E_from_r.for'
      !include 'include/pd_pur_R2/pd_pur.2.5.for'  !Included pd_dis.1.4.for
      SbEk=0.d0
      SbEp=0.d0
      SbE2 = 0.d0
      do ib = 1, nb
       !Divisibility of block index
       !include 'include/pd_pur_R2/pd_pur.3.5.for'
       SsEp = 0.d0
       SsEk = 0.d0
       ob = ib .GT. nbSkip
       ibEff = ib - nbSkip
       do is = 1, ns
        SwEp = 0.d0
        SwEk = 0.d0
        jw  =  0
        do iw = 1, nw
          Eold = Elocal(iw,iso)

************************************************************************
* 2 **    Gausssian displacement: r1=rp+x; x is randomly drawn    ******
******    from the 3N gaussian distribution exp(-x*x/(4*D*dt))    ******
          n=1
          do ip = 1, npA
            do k=1,3
              rn(k,ip,n)=rp(k,ip,iw,iso)+gasdev(nrg)*var(1)
            enddo
          enddo
          do ip = npA+1, npAB
            do k=1,3
              rn(k,ip,n)=rp(k,ip,iw,iso)+gasdev(nrg)*var(2)
            enddo
          enddo
          do ip = npAB+1, np
            do k=1,3
              rn(k,ip,n)=rp(k,ip,iw,iso)+gasdev(nrg)*var(3)
            enddo
          enddo

************************************************************************
* 3 **    Calculating half of drift force F1(R1)                  ******
          include 'include/DO_calc_R1_R2.for'
          include 'include/TWF_F-Ek.for'

************************************************************************
* 4 **    Auxliliary drift movement: R2=R1+0.5*D*dt*F1(R1)        ******
          n=2
          !There is no factor 0.5 because F is half of quantum force
          do ip = 1, npA
            do k=1,3
              rn(k,ip,n)=rn(k,ip,n-1)-Ck(1)*tau*F(k,ip,n-1)
            enddo
          enddo
          do ip = npA+1, npAB
            do k=1,3
              rn(k,ip,n)=rn(k,ip,n-1)-Ck(2)*tau*F(k,ip,n-1)
            enddo
          enddo
          do ip = npAB+1, np
            do k=1,3
              rn(k,ip,n)=rn(k,ip,n-1)-Ck(3)*tau*F(k,ip,n-1)
            enddo
          enddo

************************************************************************
* 5 **    Calculation of the drift force F2(R2)                   ******
          include 'include/DO_calc_R1_R2.for'
          include 'include/TWF_F-Ek.for'

************************************************************************
* 6 **    Drift movement to the middle point R2=R1+0.25Ddt(F1+F2) ******
          do ip = 1, npA
            do k=1,3
              rn(k,ip,n)=rn(k,ip,n-1)
     &                  -0.5d0*Ck(1)*tau*(F(k,ip,n-1)+F(k,ip,n))
            enddo
          enddo
          do ip = npA+1, npAB
            do k=1,3
              rn(k,ip,n)=rn(k,ip,n-1)
     &                  -0.5d0*Ck(2)*tau*(F(k,ip,n-1)+F(k,ip,n))
            enddo
          enddo
          do ip = npAB+1, np
            do k=1,3
              rn(k,ip,n)=rn(k,ip,n-1)
     &                  -0.5d0*Ck(3)*tau*(F(k,ip,n-1)+F(k,ip,n))
            enddo
          enddo

************************************************************************
* 7 **    Calculation of the local energy, drift force F(R2) etc. ******
          include 'include/DO_calc_R1_R2.for'
          include 'include/E_Ek.for'
          include 'include/E_Ep.for'
          Enew = SEp(iw) + SEk(iw)

************************************************************************
* 8 **    Final drift movement: R3=R1+D*dt*F2                     ******
          do ip = 1, npA
            do k=1,3
              rn(k,ip,n+1)=rn(k,ip,n-1)-2.d0*Ck(1)*tau*F(k,ip,n)
            enddo
          enddo
          do ip = npA+1, npAB
            do k=1,3
              rn(k,ip,n+1)=rn(k,ip,n-1)-2.d0*Ck(2)*tau*F(k,ip,n)
            enddo
          enddo
          do ip = npAB+1, np
            do k=1,3
              rn(k,ip,n+1)=rn(k,ip,n-1)-2.d0*Ck(3)*tau*F(k,ip,n)
            enddo
          enddo

************************************************************************
* 9 **    Branching weight: w=exp(-dt*(0.5*(ELn+ELo)-E))          ******
          if(nwRef.GT.1)then
            rsons = dexp(tau*(AwE-0.5d0*(Enew+Eold)))  
            nsons = rsons
            fsons = rsons-nsons
            if(fsons.GT.0.1d0)then
              nsons = int(rsons+ran1(nrg))
            else
              prob =dsqrt(dsqrt(dsqrt(fsons))) 
              o1 = ran1(nrg).LT.prob
              o2 = ran1(nrg).LT.prob
              o3 = ran1(nrg).LT.prob
              o4 = ran1(nrg).LT.prob
              o5 = ran1(nrg).LT.prob
              o6 = ran1(nrg).LT.prob
              o7 = ran1(nrg).LT.prob
              o8 = ran1(nrg).LT.prob
              o1 = (o1.AND.o2).AND.(o3.AND.o4)
              o2 = (o5.AND.o6).and.(o7.AND.o8)
              o  = (o1.AND.o2)
              if (o) nsons = nsons + 1
            endif
          else
            nsons=1
          endif
          if (nsons.GT.0 .AND. nwRef.NE.1)then
            if (nw.GT.nwMax)then
              xsons = nsons * redu
              nsaux = xsons + ran1(nrg)
              nsons = nsaux
            endif
            if(nw.LT.nwMin)then
               xsons = nsons * ampi 
               nsaux = xsons + ran1(nrg)
               nsons = nsaux
            endif
          endif
          !Acumulating estimators
          if (nsons.GT.0)then
            do isons = 1, nsons
              jw = jw + 1
              if(jw.GT.mw)then
                write(iScr,'(A5,A20,I8,A3,I8)')TRIM(ADJUSTL(wTau)),
     &                     ' ===>>>>>>>  ERROR: ', jw, ' > ', mw
                CLOSE(iScr)
                STOP
              endif
              Elocal(jw,isn)=Enew
              do ip = 1, np
                do k = 1, 3
                  rp(k,ip,jw,isn)=rn(k,ip,n+1)
                enddo
              enddo
              !Acumulating pur estimators
              !include 'include/pd_pur_R2/pd_pur.4.5.for' !included pd_dis.3.4.for
            enddo !isons
            SwEp = SwEp + SEp(iw)*nsons
            SwEk = SwEk + SEk(iw)*nsons
          endif
        enddo !walkers
        nw  = jw
        iso = 3 - iso
        isn = 3 - isn
        AwEp = SwEp / jw
        AwEk = SwEk / jw
        AwE  = AwEp+AwEk
        SsEp = SsEp+AwEp
        SsEk = SsEk+AwEk
       enddo !is
       AsEp = SsEp/ns
       AsEk = SsEk/ns
       if(ob)then
         SbEp = SbEp+AsEp
         SbEk = SbEk+AsEk
       endif
       include 'include/E_save_DMC.for'

       !Writing average values of pure estimators if lpur equal 1
       !include 'include/pd_pur_R2/pd_pur.5.5.for' !Included pd_dis.4m.4.for 
                                                  !Included pd_dis.4p.4.for
      enddo !ib
      AbE=(SbEk+SbEp)/ibEff
      write(iScr,'(A24)')'iTau, dTau mK, E / mK = '
      write(iScr,'(A5,2G15.8)')TRIM(ADJUSTL(wTau)), tau, AbE
      CLOSE(iScr)
      STOP
      END


      !Gaussian distribution 
      include 'include/ran_gasdev.for'
      !NR random number generator 
      include 'include/ran2.for'

