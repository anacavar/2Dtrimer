      !Screen-file
      OPEN(iScr,FILE=TRIM(wSubDir)//'Screen.inf')
      wp11='#         Ek /mK          Ep /mK           E /'
      wp21='# ---------- ---  ---------- ---  ---------- -'
      wp12='mK      sigmaE /mK Acc  step/A  parameters'
      wp22='--  ---------- ---  --  ------  -------->>'
      write(iScr,'(A88)')TRIM(ADJUSTL(wp11))//TRIM(ADJUSTL(wp12))
      write(iScr,'(A88)')TRIM(ADJUSTL(wp21))//TRIM(ADJUSTL(wp22))
      !Reading parameters of VMC simulation
      OPEN(11,FILE=TRIM(wSubDir)//'in_VMC.ini')
      read(11,*)wtmp, nw
      read(11,*)wtmp, nb
      read(11,*)wtmp, ns
      read(11,*)wtmp, nrg
      read(11,*)wtmp, nbSkip
      read(11,*)wtmp, step
      read(11,*)wtmp, delta
      read(11,*)wtmp, loadr
      CLOSE(11)
      sTemp = step
      nbTmp = nb ! needed while changing loadr during simulation
      loadrTmp = loadr
      nbSkipTmp = nbSkip
      if(loadr .EQ. 1) then
        nRun = 2
        nb = nbSkip + 50
      else
        nRun = 1
      endif
      DF=0 !To reset F if Fi doesn't exist
      
      ALLOCATE(r1(2:np,1:np-1), r2(2:np,1:np-1), dr(2,2:np,1:np-1))
      ALLOCATE(rp(2,np,nw), rn(2,np,1), Aux(2:np,1:np-1), F(2,np,1))
      !ALLOCATE(Pold(nw), Pnew(nw))
      ALLOCATE(Ek(np,nw), SEk(nw), SEp(nw))

      ALLOCATE(ePnew(nw),ePold(nw),STAT=iALLOC) !Psi parts
      if(iALLOC.NE.0) then
        write(iScr,*)'ERROR: Allocating ePsi new, old...'
        STOP
      endif
