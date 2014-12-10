#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#     config/paraview.csh
#
# Description
#     Setup file for paraview-[3-4].x
#     Sourced from OpenFOAM-<VERSION>/etc/cshrc or from foamPV alias
#
# Note
#     The env. variables 'ParaView_DIR' and 'ParaView_MAJOR'
#     are required for building plugins
#------------------------------------------------------------------------------

# clean the PATH
set cleaned=`$WM_PROJECT_DIR/bin/foamCleanPath "$PATH" "$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/cmake- $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/paraview-"`
if ( $status == 0 ) setenv PATH $cleaned

# determine the cmake to be used
unsetenv CMAKE_HOME
foreach cmake ( cmake-2.8.12.1 cmake-2.8.8 cmake-2.8.4 cmake-2.8.3 cmake-2.8.1 )
    set cmake=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$cmake
    if ( -r $cmake ) then
        setenv CMAKE_HOME $cmake
        setenv PATH ${CMAKE_HOME}/bin:${PATH}
        break
    endif
end

#- ParaView version, automatically determine major version:
#setenv ParaView_VERSION 3.12.0
#setenv ParaView_VERSION 4.0.1
setenv ParaView_VERSION 4.1.0
setenv ParaView_MAJOR detect


# Evaluate command-line parameters for ParaView
while ( $#argv > 0 )
    switch ($argv[1])
    case ParaView*=*:
        # name=value  -> setenv name value
        eval "setenv $argv[1]:s/=/ /"
        breaksw
    endsw
    shift
end


# set MAJOR version to correspond to VERSION
# ParaView_MAJOR is "<digits>.<digits>" from ParaView_VERSION
switch ("$ParaView_VERSION")
case "$ParaView_MAJOR".*:
    # version and major appear to correspond
    breaksw

case [0-9]*:
    # extract major from the version
    setenv ParaView_MAJOR `echo ${ParaView_VERSION} | \
        sed -e 's/^\([0-9][0-9]*\.[0-9][0-9]*\).*$/\1/'`
    breaksw
endsw


set paraviewInstDir=$WM_THIRD_PARTY_DIR/ParaView-${ParaView_VERSION}
set paraviewArchName=ParaView-$ParaView_VERSION

# Reset the name of the binary install directory for version 3
if ( `echo $ParaView_VERSION | sed -e 's/^\([0-9][0-9]*\).*$/\1/'` == 3) then
    set paraviewArchName=paraview-$ParaView_VERSION
endif

setenv ParaView_DIR $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$paraviewArchName

# set paths if binaries or source are present
if ( -r $ParaView_DIR || -r $paraviewInstDir ) then
    setenv ParaView_INCLUDE_DIR $ParaView_DIR/include/paraview-${ParaView_MAJOR}
    if (! -r $ParaView_INCLUDE_DIR && -r $ParaView_DIR/include/paraview) then
        setenv ParaView_INCLUDE_DIR $ParaView_DIR/include/paraview
    endif

    set ParaView_LIB_DIR=${ParaView_DIR}/lib/paraview-${ParaView_MAJOR}
    if (! -r $ParaView_LIB_DIR && -r ${ParaView_DIR}/lib/paraview) then
        set ParaView_LIB_DIR=${ParaView_DIR}/lib/paraview
    endif

    setenv PATH ${ParaView_DIR}/bin:${PATH}
    setenv LD_LIBRARY_PATH "${ParaView_LIB_DIR}:${LD_LIBRARY_PATH}"
    setenv PV_PLUGIN_PATH $FOAM_LIBBIN/paraview-${ParaView_MAJOR}

    if ($?FOAM_VERBOSE && $?prompt) then
        echo "Using paraview"
        echo "    ParaView_DIR         : $ParaView_DIR"
        echo "    ParaView_LIB_DIR     : $ParaView_LIB_DIR"
        echo "    ParaView_INCLUDE_DIR : $ParaView_INCLUDE_DIR"
        echo "    PV_PLUGIN_PATH       : $PV_PLUGIN_PATH"
    endif


    # add in python libraries if required
    set paraviewPython=$ParaView_DIR/Utilities/VTKPythonWrapping
    if ( -r $paraviewPython ) then
        if ($?PYTHONPATH) then
            setenv PYTHONPATH ${PYTHONPATH}:${paraviewPython}:$ParaView_LIB_DIR
        else
            setenv PYTHONPATH ${paraviewPython}:$ParaView_LIB_DIR
        endif
    endif
else
    unsetenv PV_PLUGIN_PATH
endif


unset cleaned cmake paraviewInstDir paraviewPython

# -----------------------------------------------------------------------------
