      !Finding system directories (wdir)
      write(wmA,*)mA
      write(wmB,*)mB
      write(wmC,*)mC
      wmA=ADJUSTL(wmA)
      wmB=ADJUSTL(wmB)
      wmC=ADJUSTL(wmC)
      !Setting work directories
      if(iDir.NE.0) then
        if(npA.GT.0) then
          wDir=TRIM(wmA)
          do i = 2, npA
            wDir=TRIM(wDir)//'.'//TRIM(wmA)
          enddo
          do i = 1, npB
            wDir=TRIM(wDir)//'.'//TRIM(wmB)
          enddo
        elseif(npB.GT.0) then
          wDir=TRIM(wmB)
          do i = 2, npB
            wDir=TRIM(wDir)//'.'//TRIM(wmB)
          enddo
        endif
        do i = 1, npC
          wDir=TRIM(wDir)//'.'//TRIM(wmC)
        enddo
        wDir=TRIM(wDir)//'/'
        wSubDir=TRIM(wDir)//TRIM(wFolder)//'/'
      else
        wDir=''
        wSubDir=''
      endif


