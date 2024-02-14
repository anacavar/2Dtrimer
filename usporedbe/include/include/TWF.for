      !Reading parameters of trial wave function (wTwf)
      ALLOCATE(a(6),g(6),s(6),STAT=iALLOC)
      if(iALLOC.NE.0) then
        write(iScr,*)'ERROR: Allocating Twf...'
        STOP
      endif
      do i = 1, 6
        a(i)=0.d0
        g(i)=0.d0
        s(i)=0.d0
      enddo
      OPEN(11,FILE=TRIM(wSubDir)//'in_TWF_ags.ini')
      read(11,*)
      if(npA.GT.1) read(11,*)a(1), g(1), s(1)
      if(npB.GT.1) read(11,*)a(2), g(2), s(2)
      if(npC.GT.1) read(11,*)a(3), g(3), s(3)
      if((npA*npB).GT.0) read(11,*)a(4), g(4), s(4)
      if((npA*npC).GT.0) read(11,*)a(5), g(5), s(5)
      if((npB*npC).GT.0) read(11,*)a(6), g(6), s(6)
      CLOSE(11)
