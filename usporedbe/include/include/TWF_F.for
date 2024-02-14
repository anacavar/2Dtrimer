! Aux-line can be commented in VMC


          !Reseting values for half of quantum force, F
          do ip = 1, np
            do k = 1, 2
              F(k,ip,n) = 0.d0
            enddo
          enddo !particles
          ! F and contribution of all pairs to Ek
          do ip = 1, np-1
            do jp = ip+1, np
              i = ijp(jp,ip)
              Aux(jp,ip) = (a(i)/r1(jp,ip))**g(i)
              pdlr = (g(i)*Aux(jp,ip)-s(i)*r1(jp,ip)-0.5d0)/r2(jp,ip)
              pddr = -1.d0*(g(i)*g(i)*Aux(jp,ip) + s(i)*r1(jp,ip))
     &             / r2(jp,ip)
              do k = 1, 2
                F(k,ip,n) = F(k,ip,n) - pdlr * dr(k,jp,ip)
                F(k,jp,n) = F(k,jp,n) + pdlr * dr(k,jp,ip)
              enddo
              Ek(ip,iw)=Ek(ip,iw)+pddr
              Ek(jp,iw)=Ek(jp,iw)+pddr
            enddo
          enddo
