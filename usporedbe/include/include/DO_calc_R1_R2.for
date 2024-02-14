          !Claculating distance beetween particles 
          do ip = 1, np-1
            do jp = ip+1, np
              dr(1,jp,ip) = rn(1,jp,n) - rn(1,ip,n)
              dr(2,jp,ip) = rn(2,jp,n) - rn(2,ip,n)
              r2(jp,ip) = dr(1,jp,ip)*dr(1,jp,ip)
     &                  + dr(2,jp,ip)*dr(2,jp,ip)
              r1(jp,ip) = dsqrt(r2(jp,ip))
            enddo
          enddo
