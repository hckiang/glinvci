module test
  use, intrinsic :: iso_c_binding
  implicit integer(c_int) (i-k), integer(c_int) (m,n), &
       & real(c_double) (a-h), real(c_double) (l), real(c_double) (o-z)
contains
  subroutine glinvtestfloatIEEE01(x,out) bind(C, name='glinvtestfloatIEEE01_')
    out=exp(x)
  end subroutine
end module
