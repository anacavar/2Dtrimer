          !Writing energies after block
          if(ibEff.EQ.1) then
            OPEN(20,FILE=TRIM(wsubdir)//'VMC_AbE.dat')
            OPEN(30,FILE=TRIM(wsubdir)//'VMC_AE.dat')
            write(20,'(A68)')'# (ib-nbSkip)*ns           E /mK'//
     &                       '          Ek /mK          Ep /mK Acc'
            write(20,'(A68)')'# --------------  ---------- ---'//
     &                       '  ---------- ---  ---------- ---  --'
            write(30,'(A68)')'# (ib-nbSkip)*ns           E /mK'//
     &                       '          Ek /mK          Ep /mK Acc'
            write(30,'(A68)')'# --------------  ---------- ---'//
     &                       '  ---------- ---  ---------- ---  --'
            CLOSE(20)
            CLOSE(30)
          endif
          OPEN(20,FILE=TRIM(wsubdir)//'VMC_AbE.dat', ACCESS='APPEND')
          OPEN(30,FILE=TRIM(wsubdir)//'VMC_AE.dat',  ACCESS='APPEND')
          write(20,'(I16,3G16.8,I4)') ibEff*ns, AbEk+AbEp, AbEk, AbEp,
     &                                NINT(Acc*100.0)
          write(30,'(I16,3G16.8,I4)') ibEff*ns, 
     &                               (SbEk+SbEp) / ibEff ,
     &                                SbEk       / ibEff ,
     &                                SbEp       / ibEff ,
     &                                NINT(Acc*100.0)
          CLOSE(20)
          CLOSE(30)
