#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     etc/config/settings.sh
#
# Description
#     Startup file for OpenFOAM
#     Sourced from OpenFOAM-<VERSION>/etc/bashrc
#
#------------------------------------------------------------------------------

# Prefix to PATH
_foamAddPath()
{
    while [ $# -ge 1 ]
    do
        export PATH=$1:$PATH
        shift
    done
}

# Prefix to LD_LIBRARY_PATH
_foamAddLib()
{
    while [ $# -ge 1 ]
    do
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
        shift
    done
}

# Prefix to MANPATH
_foamAddMan()
{
    while [ $# -ge 1 ]
    do
        export MANPATH=$1:$MANPATH
        shift
    done
}

#------------------------------------------------------------------------------
# Set environment variables according to system type
export WM_ARCH=`uname -s`

case "$WM_ARCH" in
Linux)
    WM_ARCH=linux

    # compiler specifics
    case `uname -m` in
        i686)
            export WM_ARCH_OPTION=32
        ;;

    x86_64)
        case "$WM_ARCH_OPTION" in
        32)
            export WM_COMPILER_ARCH=64
            export WM_CC='gcc'
            export WM_CXX='g++'
            export WM_CFLAGS='-m32 -fPIC'
            export WM_CXXFLAGS='-m32 -fPIC'
            export WM_LDFLAGS='-m32'
            ;;
        64)
            WM_ARCH=linux64
            export WM_COMPILER_LIB_ARCH=64
            export WM_CC='gcc'
            export WM_CXX='g++'
            export WM_CFLAGS='-m64 -fPIC'
            export WM_CXXFLAGS='-m64 -fPIC'
            export WM_LDFLAGS='-m64'
            ;;
        *)
            echo "Unknown WM_ARCH_OPTION '$WM_ARCH_OPTION', should be 32 or 64"\
                 1>&2
            ;;
        esac
        ;;

    ia64)
        WM_ARCH=linuxIA64
        export WM_COMPILER=I64
        ;;

    armv7l)
        WM_ARCH=linuxARM7
        export WM_COMPILER_LIB_ARCH=32
        export WM_CC='gcc'
        export WM_CXX='g++'
        export WM_CFLAGS='-fPIC'
        export WM_CXXFLAGS='-fPIC'
        export WM_LDFLAGS=
        ;;

    ppc64)
        WM_ARCH=linuxPPC64
        export WM_COMPILER_LIB_ARCH=64
        export WM_CC='gcc'
        export WM_CXX='g++'
        export WM_CFLAGS='-m64 -fPIC'
        export WM_CXXFLAGS='-m64 -fPIC'
        export WM_LDFLAGS='-m64'
        ;;

    ppc64le)
        WM_ARCH=linuxPPC64le
        export WM_COMPILER_LIB_ARCH=64
        export WM_CC='gcc'
        export WM_CXX='g++'
        export WM_CFLAGS='-m64 -fPIC'
        export WM_CXXFLAGS='-m64 -fPIC'
        export WM_LDFLAGS='-m64'
        ;;

    *)
        echo Unknown processor type `uname -m` for Linux 1>&2
        ;;
    esac
    ;;

SunOS)
    WM_ARCH=SunOS64
    WM_MPLIB=FJMPI
    export WM_COMPILER_LIB_ARCH=64
    export WM_CC='gcc'
    export WM_CXX='g++'
    export WM_CFLAGS='-mabi=64 -fPIC'
    export WM_CXXFLAGS='-mabi=64 -fPIC'
    export WM_LDFLAGS='-mabi=64 -G0'
    ;;

*)    # an unsupported operating system
    /bin/cat <<USAGE 1>&2

    Your "$WM_ARCH" operating system is not supported by this release
    of OpenFOAM. For further assistance, please contact www.OpenFOAM.org

USAGE
    ;;
esac


#------------------------------------------------------------------------------

# Location of the jobControl directory
export FOAM_JOB_DIR=$WM_PROJECT_INST_DIR/jobControl

# wmake configuration
export WM_DIR=$WM_PROJECT_DIR/wmake
export WM_LINK_LANGUAGE=c++
export WM_LABEL_OPTION=Int$WM_LABEL_SIZE
export WM_OPTIONS=$WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_LABEL_OPTION$WM_COMPILE_OPTION

# Base executables/libraries
export FOAM_APPBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/bin
export FOAM_LIBBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib

# External (ThirdParty) libraries
export FOAM_EXT_LIBBIN=$WM_THIRD_PARTY_DIR/platforms/$WM_OPTIONS/lib

# Site-specific directory
siteDir="${WM_PROJECT_SITE:-$WM_PROJECT_INST_DIR/site}"

# Shared site executables/libraries
# Similar naming convention as ~OpenFOAM expansion
export FOAM_SITE_APPBIN=$siteDir/$WM_PROJECT_VERSION/platforms/$WM_OPTIONS/bin
export FOAM_SITE_LIBBIN=$siteDir/$WM_PROJECT_VERSION/platforms/$WM_OPTIONS/lib

# User executables/libraries
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib

# DynamicCode templates
# - default location is the "~OpenFOAM/codeTemplates/dynamicCode" expansion
# export FOAM_CODE_TEMPLATES=$WM_PROJECT_DIR/etc/codeTemplates/dynamicCode

# Convenience
export FOAM_ETC=$WM_PROJECT_DIR/etc
export FOAM_APP=$WM_PROJECT_DIR/applications
export FOAM_SRC=$WM_PROJECT_DIR/src
export FOAM_TUTORIALS=$WM_PROJECT_DIR/tutorials
export FOAM_UTILITIES=$FOAM_APP/utilities
export FOAM_SOLVERS=$FOAM_APP/solvers
export FOAM_RUN=$WM_PROJECT_USER_DIR/run

# Add wmake to the path - not required for runtime-only environment
[ -d "$WM_DIR" ] && PATH=$WM_DIR:$PATH
# Add OpenFOAM scripts to the path
export PATH=$WM_PROJECT_DIR/bin:$PATH

# add site-specific scripts to path - only if they exist
if [ -d "$siteDir/bin" ]                        # generic
then
    _foamAddPath "$siteDir/bin"
fi
if [ -d "$siteDir/$WM_PROJECT_VERSION/bin" ]    # version-specific
then
    _foamAddPath "$siteDir/$WM_PROJECT_VERSION/bin"
fi
unset siteDir

_foamAddPath $FOAM_USER_APPBIN:$FOAM_SITE_APPBIN:$FOAM_APPBIN
# Make sure to pick up dummy versions of external libraries last
_foamAddLib  $FOAM_USER_LIBBIN:$FOAM_SITE_LIBBIN:$FOAM_LIBBIN:$FOAM_EXT_LIBBIN:$FOAM_LIBBIN/dummy

# Compiler settings
# ~~~~~~~~~~~~~~~~~
unset gcc_version gmp_version mpfr_version mpc_version
unset MPFR_ARCH_PATH GMP_ARCH_PATH

# Location of compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if [ -z "$foamCompiler" ]
then
    foamCompiler=system
    echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
    echo "    foamCompiler not set, using '$foamCompiler'" 1>&2
fi

case "${foamCompiler}" in
OpenFOAM | ThirdParty)
    case "$WM_COMPILER" in
    Gcc | Gcc48)
        gcc_version=gcc-4.8.4
        gmp_version=gmp-5.1.2
        mpfr_version=mpfr-3.1.2
        mpc_version=mpc-1.0.1
        ;;
    Gcc45)
        gcc_version=gcc-4.5.4
        gmp_version=gmp-5.1.2
        mpfr_version=mpfr-3.1.2
        mpc_version=mpc-1.0.1
        ;;
    Gcc46)
        gcc_version=gcc-4.6.4
        gmp_version=gmp-5.1.2
        mpfr_version=mpfr-3.1.2
        mpc_version=mpc-1.0.1
        ;;
    Gcc47)
        gcc_version=gcc-4.7.4
        gmp_version=gmp-5.1.2
        mpfr_version=mpfr-3.1.2
        mpc_version=mpc-1.0.1
        ;;
    Gcc49)
        gcc_version=gcc-4.9.2
        gmp_version=gmp-5.1.2
        mpfr_version=mpfr-3.1.2
        mpc_version=mpc-1.0.1
        ;;
    Gcc51)
        gcc_version=gcc-5.1.0
        gmp_version=gmp-5.1.2
        mpfr_version=mpfr-3.1.2
        mpc_version=mpc-1.0.1
        ;;
    Clang)
        # using clang - not gcc
        export WM_CC='clang'
        export WM_CXX='clang++'
        clang_version=llvm-3.6.0
        ;;
    *)
        echo 1>&2
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    Unknown OpenFOAM compiler type '$WM_COMPILER'" 1>&2
        echo "    Please check your settings" 1>&2
        echo 1>&2
        ;;
    esac

    # Optional configuration tweaks:
    _foamSource `$WM_PROJECT_DIR/bin/foamEtcFile config/compiler.sh`

    if [ -n "$gcc_version" ]
    then
        gccDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gcc_version
        gmpDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gmp_version
        mpfrDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpfr_version
        mpcDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpc_version

        # Check that the compiler directory can be found
        [ -d "$gccDir" ] || {
            echo 1>&2
            echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
            echo "    Cannot find $gccDir installation." 1>&2
            echo "    Please install this compiler version or if you wish to" \
                 " use the system compiler," 1>&2
            echo "    change the 'foamCompiler' setting to 'system'" 1>&2
            echo
        }

        _foamAddMan     $gccDir/man
        _foamAddPath    $gccDir/bin

        # Add compiler libraries to run-time environment
        _foamAddLib     $gccDir/lib$WM_COMPILER_LIB_ARCH

        # Add gmp/mpfr libraries to run-time environment
        _foamAddLib     $gmpDir/lib
        _foamAddLib     $mpfrDir/lib

        # Add mpc libraries (not need for older gcc) to run-time environment
        if [ -n "$mpc_version" ]
        then
            _foamAddLib     $mpcDir/lib
        fi

        # Used by boost/CGAL:
        export MPFR_ARCH_PATH=$mpfrDir
        export GMP_ARCH_PATH=$gmpDir
    fi
    unset gcc_version gccDir
    unset gmp_version gmpDir  mpfr_version mpfrDir  mpc_version mpcDir

    if [ -n "$clang_version" ]
    then
        clangDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$clang_version

        # Check that the compiler directory can be found
        [ -d "$clangDir" ] || {
            echo 1>&2
            echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
            echo "    Cannot find $clangDir installation." 1>&2
            echo "    Please install this compiler version or if you wish to" \
                 " use the system compiler," 1>&2
            echo "    change the 'foamCompiler' setting to 'system'" 1>&2
            echo 1>&2
        }

        _foamAddMan     $clangDir/share/man
        _foamAddPath    $clangDir/bin
    fi
    unset clang_version clangDir
    ;;
system)
    # Use system compiler
    ;;
*)
    echo "Warn: foamCompiler='$foamCompiler' is unsupported" 1>&2
    echo "   treating as 'system' instead" 1>&2
    ;;
esac


#
# Add c++0x flags for external programs
#
if [ -n "$WM_CXXFLAGS" ]
then
    case "$WM_COMPILER" in
    Gcc*++0x)
        WM_CXXFLAGS="$WM_CXXFLAGS -std=c++0x"
        ;;
    esac
fi



# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH MPI_HOME FOAM_MPI_LIBBIN

case "$WM_MPLIB" in
SYSTEMOPENMPI)
    # Use the system installed openmpi, get library directory via mpicc
    export FOAM_MPI=openmpi-system

    libDir=`mpicc --showme:link | sed -e 's/.*-L\([^ ]*\).*/\1/'`

    # Bit of a hack: strip off 'lib' and hope this is the path to openmpi
    # include files and libraries.
    export MPI_ARCH_PATH="${libDir%/*}"

    _foamAddLib     $libDir
    unset libDir
    ;;

OPENMPI)
    export FOAM_MPI=openmpi-1.10.0
    # Optional configuration tweaks:
    _foamSource `$WM_PROJECT_DIR/bin/foamEtcFile config/openmpi.sh`

    export MPI_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$FOAM_MPI

    # Tell OpenMPI where to find its install directory
    export OPAL_PREFIX=$MPI_ARCH_PATH

    _foamAddPath    $MPI_ARCH_PATH/bin

    # 64-bit on OpenSuSE 12.1 uses lib64 others use lib
    _foamAddLib     $MPI_ARCH_PATH/lib$WM_COMPILER_LIB_ARCH
    _foamAddLib     $MPI_ARCH_PATH/lib

    _foamAddMan     $MPI_ARCH_PATH/share/man
    ;;

SYSTEMMPI)
    export FOAM_MPI=mpi-system

    if [ -z "$MPI_ROOT" ]
    then
        echo 1>&2
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    Please set the environment variable MPI_ROOT to point to" \
             " the base folder for the system MPI in use." 1>&2
        echo "    Example:" 1>&2
        echo 1>&2
        echo "        export MPI_ROOT=/opt/mpi" 1>&2
        echo 1>&2
    else
        export MPI_ARCH_PATH=$MPI_ROOT

        if [ -z "$MPI_ARCH_FLAGS" ]
        then
            echo 1>&2
            echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
            echo "    MPI_ARCH_FLAGS is not set. Example:" 1>&2
            echo 1>&2
            echo "        export MPI_ARCH_FLAGS=\"-DOMPI_SKIP_MPICXX\"" 1>&2
            echo 1>&2
        fi

        if [ -z "$MPI_ARCH_INC" ]
        then
            echo 1>&2
            echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
            echo "    MPI_ARCH_INC is not set. Example:" 1>&2
            echo 1>&2
            echo "        export MPI_ARCH_INC=\"-isystem \$MPI_ROOT/include\"" 1>&2
            echo 1>&2
        fi

        if [ -z "$MPI_ARCH_LIBS" ]
        then
            echo 1>&2
            echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
            echo "    MPI_ARCH_LIBS is not set. Example:" 1>&2
            echo 1>&2
            echo "        export MPI_ARCH_LIBS=\"-L\$MPI_ROOT/lib -lmpi\"" 1>&2
            echo 1>&2
        fi
    fi
    ;;

MPICH)
    export FOAM_MPI=mpich2-1.1.1p1
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$FOAM_MPI
    export MPI_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$FOAM_MPI

    _foamAddPath    $MPI_ARCH_PATH/bin

    # 64-bit on OpenSuSE 12.1 uses lib64 others use lib
    _foamAddLib     $MPI_ARCH_PATH/lib$WM_COMPILER_LIB_ARCH
    _foamAddLib     $MPI_ARCH_PATH/lib

    _foamAddMan     $MPI_ARCH_PATH/share/man
    ;;

MPICH-GM)
    export FOAM_MPI=mpich-gm
    export MPI_ARCH_PATH=/opt/mpi
    export MPICH_PATH=$MPI_ARCH_PATH
    export GM_LIB_PATH=/opt/gm/lib64

    _foamAddPath    $MPI_ARCH_PATH/bin

    # 64-bit on OpenSuSE 12.1 uses lib64 others use lib
    _foamAddLib     $MPI_ARCH_PATH/lib$WM_COMPILER_LIB_ARCH
    _foamAddLib     $MPI_ARCH_PATH/lib

    _foamAddLib     $GM_LIB_PATH
    ;;

HPMPI)
    export FOAM_MPI=hpmpi
    export MPI_HOME=/opt/hpmpi
    export MPI_ARCH_PATH=$MPI_HOME

    _foamAddPath $MPI_ARCH_PATH/bin

    case `uname -m` in
    i686)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia32
        ;;

    x86_64)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_amd64
        ;;
    ia64)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia64
        ;;
    *)
        echo Unknown processor type `uname -m` 1>&2
        ;;
    esac
    ;;

MPI)
    export FOAM_MPI=mpi
    export MPI_ARCH_PATH=/opt/mpi
    ;;

FJMPI)
    export FOAM_MPI=fjmpi
    export MPI_ARCH_PATH=/opt/FJSVmpi2

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib/sparcv9
    _foamAddLib     /opt/FSUNf90/lib/sparcv9
    _foamAddLib     /opt/FJSVpnidt/lib
    ;;

QSMPI)
    export FOAM_MPI=qsmpi
    export MPI_ARCH_PATH=/usr/lib/mpi

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    ;;

SGIMPI)
    # No trailing slash
    [ "${MPI_ROOT%/}" = "${MPI_ROOT}" ] || MPI_ROOT="${MPI_ROOT%/}"

    export FOAM_MPI="${MPI_ROOT##*/}"
    export MPI_ARCH_PATH=$MPI_ROOT

    if [ ! -d "$MPI_ROOT" -o -z "$MPI_ARCH_PATH" ]
    then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    MPI_ROOT not a valid mpt installation directory or ending" \
             " in a '/'." 1>&2
        echo "    Please set MPI_ROOT to the mpt installation directory." 1>&2
        echo "    MPI_ROOT currently set to '$MPI_ROOT'" 1>&2
    fi

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using SGI MPT:" 1>&2
        echo "    MPI_ROOT : $MPI_ROOT" 1>&2
        echo "    FOAM_MPI : $FOAM_MPI" 1>&2
    fi

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    ;;

INTELMPI)
    # No trailing slash
    [ "${MPI_ROOT%/}" = "${MPI_ROOT}" ] || MPI_ROOT="${MPI_ROOT%/}"

    export FOAM_MPI="${MPI_ROOT##*/}"
    export MPI_ARCH_PATH=$MPI_ROOT

    if [ ! -d "$MPI_ROOT" -o -z "$MPI_ARCH_PATH" ]
    then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    MPI_ROOT not a valid mpt installation directory or ending" \
             " in a '/'." 1>&2
        echo "    Please set MPI_ROOT to the mpt installation directory." 1>&2
        echo "    MPI_ROOT currently set to '$MPI_ROOT'" 1>&2
    fi

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using INTEL MPI:" 1>&2
        echo "    MPI_ROOT : $MPI_ROOT" 1>&2
        echo "    FOAM_MPI : $FOAM_MPI" 1>&2
    fi

    _foamAddPath    $MPI_ARCH_PATH/bin64
    _foamAddLib     $MPI_ARCH_PATH/lib64
    ;;
*)
    export FOAM_MPI=dummy
    ;;
esac

# Add (non-dummy) MPI implementation
# Dummy MPI already added to LD_LIBRARY_PATH and has no external libraries
if [ "$FOAM_MPI" != dummy ]
then
    _foamAddLib $FOAM_LIBBIN/$FOAM_MPI:$FOAM_EXT_LIBBIN/$FOAM_MPI
fi



# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
: ${minBufferSize:=20000000}

if [ "${MPI_BUFFER_SIZE:=$minBufferSize}" -lt $minBufferSize ]
then
    MPI_BUFFER_SIZE=$minBufferSize
fi
export MPI_BUFFER_SIZE


# Cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
#keep _foamAddPath _foamAddLib _foamAddMan
unset foamCompiler minBufferSize

# ----------------------------------------------------------------- end-of-file
