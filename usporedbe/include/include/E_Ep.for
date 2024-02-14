          !Potential energy
          SEp(iw)=0.d0
          ! Potential energy
          do ip = 1, np-1
            do jp = ip+1, np
              i = ijp(jp,ip)
              sr = sig(i) * sig(i) / r2(jp,ip)
              sr6 = sr * sr * sr
              SEp(iw) = SEp(iw) + eps4(i) * sr6 * (sr6 - 1.d0) !mK
            enddo
          enddo
