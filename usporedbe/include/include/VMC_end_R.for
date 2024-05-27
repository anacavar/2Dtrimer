      !Writing final position ...
      if(loadr.EQ.1) then
        wrep = 'VMC_R.ini'
      else
        wrep = 'VMC_R.dat'
      endif
      OPEN (11,FILE='./'//TRIM(wSubDir)//TRIM(wrep))
      write(11,'(A10,2I5)')'# nw; np: ', nw, np
      do iw = 1, nw
        do ip = 1, np
          write(11,*)rp(1,ip,iw), rp(2,ip,iw)
        enddo
        write(11,*)
        write(11,*)
      enddo
      CLOSE(11)
