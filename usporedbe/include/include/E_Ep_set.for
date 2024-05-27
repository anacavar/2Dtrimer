          do i = 1, 6
            eps4(i) = 0.d0
            sig(i)  = 0.d0
          enddo

          OPEN(11,FILE=TRIM(wSubDir)//'in_Ep.ini')
          read(11,*)
          if(npA.GT.1) read(11,*)sig(1), eps4(1)
          if(npB.GT.1) read(11,*)sig(2), eps4(2)
          if(npC.GT.1) read(11,*)sig(3), eps4(3)
          if((npA*npB).GT.0) read(11,*)sig(4), eps4(4)
          if((npA*npC).GT.0) read(11,*)sig(5), eps4(5)
          if((npB*npC).GT.0) read(11,*)sig(6), eps4(6)
          CLOSE(11)
          do ip = 1, 6
            eps4(ip) = 4.d0 * eps4(ip)
          enddo

          !Pair types
          ALLOCATE(ijp(2:np,1:np-1),STAT=iALLOC)
          if(iALLOC.NE.0) then
            write(iScr,*)'ERROR: Allocating pair types...'
            STOP
          endif
          do ip = 1, npA-1 ! AA pairs
            do jp = ip+1, npA
              ijp(jp,ip)=1
            enddo
          enddo
          do ip = npA+1, npAB-1 ! BB pairs
            do jp = ip+1, npAB
              ijp(jp,ip)=2
            enddo
          enddo
          do ip = npAB+1, np-1 ! CC pairs
            do jp = ip+1, np
              ijp(jp,ip)=3
            enddo
          enddo
          do ip = 1, npA ! AB pairs
            do jp = npA+1, npAB
              ijp(jp,ip)=4
            enddo
          enddo
          do ip = 1, npA ! AC pairs
            do jp = npAB+1, np
              ijp(jp,ip)=5
            enddo
          enddo
          do ip = npA+1, npAB ! BC pairs
            do jp = npAB+1, np
              ijp(jp,ip)=6
            enddo
          enddo

