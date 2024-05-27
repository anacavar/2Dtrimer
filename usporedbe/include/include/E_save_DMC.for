          !Writing energies after block
          if(ibEff.EQ.1) then
            OPEN(11,FILE=TRIM(wSubDir)//'DMC_AbE_'//TRIM(wTau)//'.dat')
            write(11,'(A74)')
     &              '# ibEff*ns    Average  E /mK    Average Ek /mK'//
     &                        '    Average Ep /mK   walkers'
            write(11,'(A74)')
     &              '# --------  ------------ ---  ------------ ---'//
     &                        '  ------------ ---  --------'
            CLOSE(11)
          endif
          OPEN(11,FILE=TRIM(wSubDir)//'DMC_AbE_'//TRIM(wTau)//'.dat',
     &            ACCESS='APPEND')
          write(11,'(I10,3G18.9,I10)') ibEff*ns,AsEk+AsEp,AsEk,AsEp,nw
          CLOSE(11)
          OPEN(11,FILE=TRIM(wSubDir)//'DMC_R_'//TRIM(wTau)//'.dat')
          write(11,*)'# nw  np; nw*(np*(x,y,z) + \n\n)'
          write(11,'(A2,2I6)')'# ', nw, np
          do iw=1,nw
            do ip=1,np
             write(11,*)rp(1,ip,iw,iso),rp(2,ip,iw,iso)
            enddo
            write(11,*)
            write(11,*)
          enddo
          CLOSE(11)
