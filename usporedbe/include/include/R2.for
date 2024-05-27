          !Potential energy3
          do i = 1, 6
            SR2(i,iw)=0.d0
          enddo
          ! Potential energy
          do ip = 1, np-1
            do jp = ip+1, np
              i = ijp(jp,ip)
              SR2(i,iw) = SR2(i,iw) + r2(jp,ip) !A**2
            enddo
          enddo
