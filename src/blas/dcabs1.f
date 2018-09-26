      double precision function dcabs1(z)
      double complex z
      double precision t(2)
      t=transfer(z,t)
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
