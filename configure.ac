AC_PREREQ(2.61)

AC_INIT(glinvci, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))
AC_COPYRIGHT(Copyright (C) 2021 Hao Chi Kiang)


## Determine Install Location of R
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi

## Get CC from R, and use this to get OpenMP flag from Autoconf.
## Note that it is impossible to get SHLIB_OPENMP_CFLAGS from R; nor is it possible
## to rely on R CMD SHLIB to supply the right OpenMP flags. So during the
## configuration stage this is the best we can do...
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
AC_LANG(C)
AC_PROG_CC([${CC}])
AC_OPENMP

## Now get some informations from R for compiling our test program
CC=`"${RBIN}" CMD config CC`
FC=`"${RBIN}" CMD config FC`
CFLAGS=`"${RBIN}" CMD config CFLAGS`
FCFLAGS=`"${RBIN}" CMD config FCFLAGS`
RBLAS=`"${RBIN}" CMD config BLAS_LIBS`
HACKY_OMPCFLAGS="${OPENMP_CFLAGS}"
RFLIBSFLAGS=`"${RBIN}" CMD config FLIBS`
RSAFEFFLAGS=`"${RBIN}" CMD config SAFE_FFLAGS`

## If autoconf has found a working OpenMP flag than we compile
## a .c file and a .f90 file, and link them together. The C
## file has OpenMP threading which calls Fortran, which in turn
## calls BLAS.
ac_pkg_openmp=no
if test -n "${OPENMP_CFLAGS}"; then
  AC_LANG_CONFTEST([AC_LANG_SOURCE([[
    #include <omp.h>
    #include <R.h>
    #include <Rinternals.h>
    #include <stdlib.h>
    #include <string.h>
    extern void fort_blastest_(int*, double*, double*, double*);
    static unsigned int g_seed = 12345U;
    static int siz = 250, nloops = 700, correct = 1, *ret;
    static double *Ac = 0;
    #pragma omp threadprivate(Ac)
    void myrand(double *x) {
      *x=((double)((((g_seed = (214013*g_seed+2531011))>>16)&0x7FFF)+1))/((double)32768);
    }
    SEXP rblastest_ABCDE() {
      double *A, *B, *Bc, *Rc, *res, fac, rsum=0.0;
      int k, j;
      Rc=0;
      omp_set_num_threads(12);
      B=(A=malloc(2*siz*siz*sizeof(double)))+siz*siz;
      for (j=0; j<siz*siz; ++j) {myrand(A+j); myrand(B+j);}
      fort_blastest_(&siz, A, B, (res=malloc(siz*siz*sizeof(double))));
      for (j=0; j<siz*siz; ++j) rsum += res[j]/siz/siz;
      #pragma omp parallel
      {
        memcpy((Ac=malloc(siz*siz*sizeof(double))),A,siz*siz*sizeof(double));
        #pragma omp barrier
        #pragma omp master
        {
          free(res);
          for (k=0; k<nloops; ++k) {
            memcpy((Bc=malloc(siz*siz*sizeof(double))),B,siz*siz*sizeof(double));
            fac = 1.0+1.5*((double)(k+1))/((double)nloops);
            for (j=0; j<siz*siz; ++j) Bc[j] *= fac;
            #pragma omp task firstprivate(Bc,fac,Rc) shared(siz,rsum,correct)
            {
              double x = (double)0.0;
              fort_blastest_(&siz, Ac, Bc, (Rc=malloc(siz*siz*sizeof(double))));
              for (int j=0; j<siz*siz; ++j) x += Rc[j]/siz/siz;
              #pragma omp critical
              correct *= (int)((x-fac*rsum)/rsum<0.00001&&(x-fac*rsum)/rsum>-0.00001);
              free(Bc); free(Rc);
            }
          }
          #pragma omp taskwait
        }
        #pragma omp barrier
        free(Ac);
      }
      SEXP Rret = PROTECT(allocVector(INTSXP, 1));
      ret = INTEGER(Rret);
      *ret = (int)(!!!correct);
      UNPROTECT(1);
      return(Rret);
    }
  ]])])
  AC_LANG_PUSH([Fortran])
  AC_FC_SRCEXT(f90)
  AC_LANG_CONFTEST([AC_LANG_SOURCE([[
module conftestfort
  use, intrinsic :: iso_c_binding
  implicit integer(c_int) (i-k), integer(c_int) (m,n), &
       & real(c_double) (a-h), real(c_double) (l), real(c_double) (o-z)
contains
  recursive subroutine fort_blastest_(k,A,B,R) bind(C,name="fort_blastest_")
    integer(c_int) :: k
    real(c_double) :: A(k,k), B(k,k), R(k,k)
    real(c_double), allocatable :: TMP(:,:)
    allocate(TMP(k,k))
    TMP(:,:) = 0.0_c_double
    call dgemm('N','N',k,k,k,1.0_c_double,A,k,B,k,0.0_c_double,TMP,k)
    call dgemm('N','N',k,k,k,1.0_c_double,A,k,TMP,k,0.0_c_double,R,k)
    deallocate(TMP)
  end subroutine
end module
  ]])])
  AC_LANG_POP([Fortran])
  ##
  ## Okay, now we've got both the C and Fortran source files; let's compile it. But
  ## 
  ## 1. We'll need to move conftest.c to some other names because
  ##      R CMD SHLIB conftest.c conftest.f90
  ##    will make R do this: gcc -shared conftest.o conftest.o.
  ##    This is an R bug.
  ##
  ## 2. Even if in the package building phase R decides to define _OPENMP and provide
  ##    an non-empty SHLIB_OPENMP_CFLAGS, the command `R CMD SHLIB' may not use -fopenmp
  ##    flag by itself. To force `R CMD SHLIB' to use -fopenmp we need to set the
  ##    PKG_CFLAGS and PKG_LIBS environment variable.
  ##
  ## 3. Windows users usually have working BLAS/OpenMP combination unless they do
  ##    something funny themselves.
  ##
  ## 4. This test has detected Rhub CentOS's BLAS being non-thread-safe! That is why
  ##    testthat fails in `R CMD check' on the platform.
  ## 
  ## 5. If we compile the test program using the $CC and $FC directly and replace the
  ##    R entry point `rblastest_ABCDE' with `int main()' it will still work
  ##    on most platforms, except Solaris. On Solaris, libRblas.so and libRlapack.so
  ##    wants a bunch of other symbols in other libraries, like -lm or -lquadmath
  ##    or some other things that I don't know of.
  ## 
  AC_MSG_CHECKING([whether BLAS works correctly with OpenMP])
  mv conftest.c conftest_c.c
  ##$CC $CFLAGS $HACKY_OMPCFLAGS -c conftest.c -o conftest_c.o && $FC $RSAFEFFLAGS $FCFLAGS -c conftest.f90 -o conftest_f.o && $CC $CFLAGS conftest_c.o conftest_f.o $HACKY_OMPCFLAGS $RLAPACK $RBLAS $RFLIBSFLAGS $RLIB -o ./conftest.exe && ./conftest.exe 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  export PKG_CFLAGS="${HACKY_OMPCFLAGS}"
  export PKG_LIBS="${RBLAS} ${HACKY_OMPCFLAGS}"
  "$RBIN" CMD SHLIB conftest_c.c conftest.f90 -o conftest.so 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep='')); quit(save='no', status=.Call('rblastest_ABCDE'), runLast=FALSE)" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  AC_MSG_RESULT([${ac_pkg_openmp}])
fi

ac_odsc=no
AC_MSG_CHECKING([whether we are on Oracle C compiler])
## Oracle Developer Studio on Solaris fails with gFortran...
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
#ifndef __SUNPRO_C
  choke me
#endif
  ])], [ac_odsc=yes])
AC_MSG_RESULT([${ac_odsc}])


if test "${ac_odsc}" = yes; then
  PRAG1="#undef __PRAGMA__"
  PRAG2="#define __PRAGMA__(x)"
  OMPFLAGS="-xopenmp=none"
  XTRA_CFLAGS="-xO1"
else
  if test "${ac_pkg_openmp}" = no; then
    if test "${OPENMP_CFLAGS}"; then
      AC_MSG_WARN([OpenMP works on your C compiler but your default BLAS obtained by "R CMD config BLAS_LIBS" has failed the thread-safety test. Turned off OpenMP multi-threading.])
    else
      AC_MSG_WARN([OpenMP is turned off manually or it does not work with your C compiler. Turning off OpenMP multi-threading.])
    fi
    PRAG1="#undef __PRAGMA__"
    PRAG2="#define __PRAGMA__(x)"
    OMPFLAGS=""
    XTRA_CFLAGS=""
  else
    PRAG1=""
    PRAG2=""
    OMPFLAGS='$(SHLIB_OPENMP_CFLAGS)'
    XTRA_CFLAGS=""
  fi  
fi

AC_SUBST(MYPRAGMA1, [${PRAG1}])
AC_SUBST(MYPRAGMA2, [${PRAG2}])
AC_SUBST(OMPFLAGS, [${OMPFLAGS}])
AC_SUBST(XTRA_CFLAGS, [${XTRA_CFLAGS}])
AC_CONFIG_FILES([src/my_pragma.h src/Makevars])
AC_OUTPUT

