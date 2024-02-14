          !Writing energies after block
          if(ibEff.EQ.1) then
            OPEN(20,FILE=TRIM(wsubdir)//'VMC_AbR2.dat')
            OPEN(30,FILE=TRIM(wsubdir)//'VMC_AR2.dat')
            write(20,'(A10,3A40)')'#    ib*ns',
     &                     '         R2(A-A)/ A2         R2(B-B)/ A2',
     &                     '         R2(C-C)/ A2         R2(A-B)/ A2',
     &                     '         R2(A-C)/ A2         R2(B-C)/ A2'
            write(20,'(A10,3A40)')'# --------',
     &                     '  -------------- ---  -------------- ---',
     &                     '  -------------- ---  -------------- ---',
     &                     '  -------------- ---  -------------- ---'
            write(30,'(A10,3A40)')'#    ib*ns',
     &                     '         R2(A-A)/ A2         R2(B-B)/ A2',
     &                     '         R2(C-C)/ A2         R2(A-B)/ A2',
     &                     '         R2(A-C)/ A2         R2(B-C)/ A2'
            write(30,'(A10,3A40)')'# --------',
     &                     '  -------------- ---  -------------- ---',
     &                     '  -------------- ---  -------------- ---',
     &                     '  -------------- ---  -------------- ---'
            CLOSE(20)
            CLOSE(30)
          endif
          OPEN(20,FILE=TRIM(wsubdir)//'VMC_AbR2.dat', ACCESS='APPEND')
          OPEN(30,FILE=TRIM(wsubdir)//'VMC_AR2.dat',  ACCESS='APPEND')
          write(20,'(I10, 6E20.12)')ibEff*ns, 
     &      AbR2(1), AbR2(2), AbR2(3), AbR2(4), AbR2(5), AbR2(6)
          write(30,'(I10, 6E20.12)')ibEff*ns, 
     &      SbR2(1)/ibEff, SbR2(2)/ibEff, SbR2(3)/ibEff, 
     &      SbR2(4)/ibEff, SbR2(5)/ibEff, SbR2(6)/ibEff
          CLOSE(20)
          CLOSE(30)
