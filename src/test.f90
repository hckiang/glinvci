module oumods
  use, intrinsic :: iso_c_binding
  implicit integer(c_int) (i-k), integer(c_int) (m,n), &
       & real(c_double) (a-h), real(c_double) (l), real(c_double) (o-z)
contains
    subroutine zI12(t,c,alpha,beta,r) bind(C, name='zI1_')
    complex(c_double_complex) c,r,    x,y,z
    if (mod2small(c) == 1) then
       r = beta*r+alpha*cmplx((t**2.0_c_double)/2.0_c_double,0.0_c_double,kind(1._c_double))
    else
       z = c*t
       y = exp(z)
       x = (2.0_c_double * cosh((z - cmplx(0._c_double,3.14159265358979324_c_double, kind(1._c_double))) &
            & /2._c_double )) / (c / exp((z + cmplx(0._c_double,3.14159265358979324_c_double,kind(1._c_double)))/2._c_double))
       r = beta*r+alpha*((t*y - x)/c)
    end if
  end subroutine
end module
