#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([OmicABELnoMM], [0.1.0], [genabel-devel@r-forge.wu-wien.ac.at])
AM_INIT_AUTOMAKE([])
AM_INIT_AUTOMAKE([silent-rules subdir-objects])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([src/Utility.h])
AC_CONFIG_HEADERS([src/config.h])

# Add the --disable-maintainer-mode option used by e.g. Debian and
# Gentoo to disable rebuilding the autotools-generated files after
# e.g. a patch has been applied to them. The default is 'disabled', we
# change it here to keep the old behaviour by default.
AM_MAINTAINER_MODE([enable])

# Set some default compile flags
if test -z "$CXXFLAGS"; then
   # User did not set CXXFLAGS, so we can put in our own defaults

   CXXFLAGS=" -O3 -std=c++11 -march=native -mfpmath=sse -flto -funroll-loops"
  #CXXFLAGS="-g -ggdb -std=c++11"

fi
if test -z "$CPPFLAGS"; then
   # User did not set CPPFLAGS, so we can put in our own defaults
   CPPFLAGS="-Wall -pedantic -Wunused-result -Wmaybe-uninitialized -Wformat"
   #CPPFLAGS="-O0"
fi
# If CXXFLAGS/CPPFLAGS are already set AC_PROG_CXX will not overwrite them
# with its own defaults


# If --with-mpi=auto is used, try to find MPI, but use standard C
#  compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi,
   [AS_HELP_STRING([--with-mpi],
      [compile with MPI (distributed memory parallelization) support. If\
 none is found, MPI is not used. Default: auto])
   ]
   ,
   ,
   [with_mpi=auto])

# AC_PROG_CXX_MPI will run AC_PROG_CXX if no MPI compiler can be found.
AX_PROG_CXX_MPI([test "x$with_mpi" != "xno"],
   [use_mpi=yes],
   [use_mpi=no
    if test "x$with_mpi" = "xyes"; then
       AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])
    elif test "x$with_mpi" = "xauto"; then
       AC_MSG_NOTICE([No MPI compiler found, continuing without MPI.])
    else
       AC_MSG_NOTICE([Use of MPI disabled by user.])
    fi]
)

AM_CONDITIONAL([USE_MPI], test "x$use_mpi" = "xyes")


# Check for openMP. If found the OPENMP_CXXFLAGS is set automatically
AC_OPENMP
AC_SUBST(AM_CXXFLAGS, "$OPENMP_CFLAGS")

#AM_CXXFLAGS="-static  -g -ggdb -I../libs/include -I./libs/include $AM_CXXFLAGS"
AM_CXXFLAGS=" -I../libs/include -I./libs/include $AM_CXXFLAGS"
# Checks for libraries.
# pthread library
AC_SEARCH_LIBS([pthread_mutex_init], [pthread], [],
   [AC_MSG_ERROR([Make sure pthread is available on the system])]
)

if test -z "$LDFLAGS"; then
   LDFLAGS="-L./libs/lib -L../libs/lib"
fi

found_blas=0

AC_SEARCH_LIBS([__iso_c_binding_c_f_pointer_l4],[gfortran])

# ACML
AC_SEARCH_LIBS(dgemm, acml_mp,
   [found_blas=1 AM_CXXFLAGS="-D_acml_ $AM_CXXFLAGS"],
   [AC_MSG_NOTICE([NOT using AMD the ACML  library])],
   [-lgfortran]
)


if test "$found_blas" -eq 0; then
   # Openblas
   AC_SEARCH_LIBS([cblas_sgemm], [openblas],
      [AM_CXXFLAGS="$AM_CXXFLAGS -D_openblas_"],
      [AC_MSG_ERROR([OpenBLAS library NOT found])]
   )
   #Lapack
   AC_SEARCH_LIBS([LAPACKE_sgeqrf], [lapack, lapacke],
      [],
      [AC_MSG_ERROR([Unable to find a Lapack library])]
   )
   AC_MSG_NOTICE([Using OpenBLAS library])
fi



# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL # USE THIS WITH AUTOCONF 2.69
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime gettimeofday memset])


# Add an option to ./configure that allows the user to specify a path
# to the Boost headers, also include the default Boost path to the
# compiler flags.
AC_ARG_WITH([boost-include-path],
   [AS_HELP_STRING([--with-boost-include-path],
      [location of the Boost math headers, defaults to /usr/include/boost/]
      )
    ],
   [CPPFLAGS+=" -I${withval}"],
   [CPPFLAGS+=' -I/usr/include/boost/']
)

# Check for the Boost Math header files; this uses a compile check, so
# make sure C++ is selected first.
AC_LANG_PUSH([C++])
AC_CHECK_HEADER([math/distributions.hpp],
  [],
  [AC_MSG_ERROR([Could not find the Boost Math header files, did you \
specify the --with-boost-include-path option correctly?])]
)
AC_LANG_POP([C++])


# The user can disable building of the PDF of the manual, for example
# if the required LaTeX packages are not installed
AC_ARG_ENABLE([latex-doc],
    [AS_HELP_STRING([--disable-latex-doc], [disable building the PDF
    documentation from LaTeX source])],
    [latexdoc=no],
    [latexdoc=yes])

if test "x$latexdoc" = "xyes"; then
   AC_MSG_NOTICE([building of the PDF of the user manual from LaTeX source is enabled])
   # Check for presence of pdfLaTeX
   AC_CHECK_PROGS(PDFLATEX, pdflatex)
   if test -z "$PDFLATEX"; then
      AC_MSG_NOTICE([pdflatex not found - Unable to create PDF version of the user manual])
   fi
fi
AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")
AM_CONDITIONAL([BUILD_latexdoc], test "x$latexdoc" = "xyes")


# Files to be generated by autotools
AC_CONFIG_FILES([
   Makefile
])

AC_OUTPUT
