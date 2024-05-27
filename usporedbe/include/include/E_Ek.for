          do ip = 1, np
            Ek(ip,iw) = 0.d0
          enddo
          include 'include/TWF_F.for' !Calculating F
          !Contribution of F*F to local energy
          do ip = 1, np
            do k = 1, 2
              Ek(ip,iw) = Ek(ip,iw) + F(k,ip,n) * F(k,ip,n)
            enddo
          enddo
          !Total Ek
          SEk(iw)=0.d0
          do ip = 1, npA
            SEk(iw) = SEk(iw) + Ek(ip,iw)*ck(1)
          enddo
          do ip = npA+1, npAB
            SEk(iw) = SEk(iw) + Ek(ip,iw)*ck(2)
          enddo
          do ip = npAB+1, np
            SEk(iw) = SEk(iw) + Ek(ip,iw)*ck(3)
          enddo
