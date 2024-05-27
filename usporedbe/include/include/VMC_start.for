      !Choosing random positions or reading them from file
      if(loadr.LT.2)then
        do iw = 1, nw
          do ip = 1, np
            rp(1,ip,iw)=delta*(ran1(nrg)-0.5d0)*2.d0
            rp(2,ip,iw)=delta*(ran1(nrg)-0.5d0)*2.d0
          enddo
        enddo
      elseif (loadr.LT.4)then  
        OPEN(11,FILE='./'//TRIM(wSubDir)//'VMC_R.ini')
        read(11,*)
        do iw = 1, nw
          do ip = 1, np
            read(11,*)rp(1,ip,iw), rp(2,ip,iw)
          enddo
          read(11,*)
          read(11,*)
        enddo
        CLOSE(11)
        !if loadr=3 do not walk (do not change position)
        if(loadr.EQ.3)then
          nbSkip=0
          step=0.d0
          ns=5 
        endif
      else
        write(iScr,*)'ERROR: loadr must be set to 0, 1, 2 or 3'
        STOP
      endif
