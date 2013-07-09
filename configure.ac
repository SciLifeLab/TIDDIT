dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(FindTranslocations, 1.0)

AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])


AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE()
LT_INIT

# trick to remove "-g -O2" default options for C*FLAGS
CXXFLAGS=""
CFLAGS=""


# standard system requirements
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_HEADERS( [pthread.h stdio.h math.h limits.h zlib.h iostream fstream string vector exception sstream], [], [AC_MSG_ERROR(["Header not found"])] )
AC_CHECK_HEADERS( [bzlib.h],[bzlib=yes] )
AC_PROG_INSTALL
AC_PROG_MAKE_SET
# AC_PROG_RANLIB

# common flags
# common flags
if test "$CXX" == "g++" ; then
	CXXFLAGS+=" -Wall -Wextra"
	CFLAGS+=" -Wall -Wextra"
fi
FLAGS+=" -D_USE_KNETFILE -D_FILE_OFFSET_BITS=64"

# checking for 64 bit system (optional)
echo $ECHO_N "checking for 64 bit architecture... "
case "$build_cpu" in
x86_64*)
	CXXFLAGS+=" -m64 -DCPU_64"
	CFLAGS+=" -m64 -DCPU_64"
	COMPILING_64="yes"
	echo "yes"
	;;
*) echo "no" ;;
esac

# checking for BOOST library (required)
echo "===== BOOST support (required) ====="
BOOST_REQUIRE([1.40.0])
BOOST_STATIC()
BOOST_FILESYSTEM
BOOST_IOSTREAMS
BOOST_PROGRAM_OPTIONS
BOOST_THREADS
BOOST_REGEX


# enable debug programs
AC_ARG_ENABLE([maintainer-mode],
     [AS_HELP_STRING([--enable-maintainer-mode],[create debug programs [no]], [28])],
     [enable_debug=$enableval],
     [enable_debug=no])
AM_CONDITIONAL(DEBUG_ENABLED,test "$enable_debug" != "no")

# enable static programs
AC_ARG_ENABLE([static-FindTranslocations],
     [AS_HELP_STRING([--enable-static-FindTranslocations],[create static FindTranslocations version [no]], [28])],
     [enable_static_FindTranslocations=$enableval],
     [enable_static_FindTranslocations=no])
AM_CONDITIONAL(STATIC_ENABLED,test "$enable_static_FindTranslocations" != "no")

AM_CONDITIONAL(GCC_COMPILER,test "$CC" = "gcc")
AM_CONDITIONAL(BGX_COMPILER,test "$CC" = "bgxlc_r")


# print summary
echo         "================ Summary ================="
if test "$COMPILING_64" = "yes" ; then
	echo "= Compiling with 64 bit                  =";
else
	echo "= Compiling with 32 bit                  =";
fi
if test "$enable_debug" != "no" ; then
	echo "= Compiling with debug programs          =";
else
	echo "= Compiling without debug programs       =";
fi
if test "$enable_static_FindTranslocations" != "no" ; then
	echo "= Compiling static FindTranslocations version  =";
else
	echo "= Not compiling static FindTranslocations version       =";
fi
echo         "=========================================="


AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT

