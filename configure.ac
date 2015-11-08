AC_INIT([gwas-pw], [0.02], [joepickrell@gmail.com], [gwas-pw], [https://github.com/joepickrell/gwas-pw])
AC_PROG_CXX
AM_INIT_AUTOMAKE([foreign tar-pax -Wall -Werror])
AC_CONFIG_HEADERS([config.h])
AC_CHECK_LIB([z],[gzopen],,AC_MSG_ERROR(could not find zlib))
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,AC_MSG_ERROR(could not find GSL BLAS))
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,AC_MSG_ERROR(could not find GSL))
m4_include([m4_ax_boost_base.m4.m4])
AX_BOOST_BASE([1.42],, AC_MSG_ERROR(need boost >= 1.42))
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS" 
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
