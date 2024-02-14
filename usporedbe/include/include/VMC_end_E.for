      !Writing final average values of energies
      if(loadr.NE.1) then
        AbEk = SbEk / ibEff
        AbEp = SbEp / ibEff
        AbE  = AbEk + AbEp
        SbAcc  = SbAcc / ibEff
        StDevE=dsqrt(dabs((SbE2/ibEff-AbE**2)/(ibEff*1.d0)))
        do i =1, 6
          write(wTwf(i),'(3F8.4)')a(i), g(i), s(i)
        enddo
        write(iScr,'(4G16.8,I4,F8.3,6A24)')
     &    AbEk, AbEp, AbE, StDevE, NINT(Acc*100.0),step,(wTwf(i),i=1,6)
      endif
