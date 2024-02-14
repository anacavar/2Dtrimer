          !Contributions
          ePnew(iw) = 0.d0
          do ip = 1, np-1
            do jp = ip+1, np
              i = ijp(jp,ip)
              Aux(jp,ip) = (a(i)/r1(jp,ip))**g(i)
              ePnew(iw) = ePnew(iw) - Aux(jp,ip) - s(i)*r1(jp,ip)
     &                              - DLOG(DSQRT(r1(jp,ip)))
            enddo
          enddo

