        SwEp = 0.d0
        SwEk = 0.d0
        n = 1
        do iw = 1, nw
          do ip = 1, np
            do k=1,2
              rn(k,ip,n)=rp(k,ip,iw,iso)
            enddo
          enddo
          include 'include/DO_calc_R1_R2.for'
          include 'include/E_Ep.for'
          include 'include/E_Ek.for'
          !Summation of energies
          Enew = SEp(iw) + SEk(iw)
          Elocal(iw,iso)=Enew
          SwEp = SwEp + SEp(iw)
          SwEk = SwEk + SEk(iw)
        enddo !walkers
        AwEp = SwEp / nw
        AwEk = SwEk / nw
        AwE  = AwEp+AwEk
