#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
#     config/CGAL.sh
#
# Description
#     Setup file for CGAL (& boost) include/libraries.
#     Sourced from OpenFOAM-<VERSION>/etc/bashrc
#------------------------------------------------------------------------------

boost_version=boost-system
cgal_version=CGAL-4.3

export BOOST_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$boost_version
export CGAL_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$cgal_version

if [ "$FOAM_VERBOSE" -a "$PS1" ]
then
    echo "Using CGAL and boost"
    echo "    $cgal_version at $CGAL_ARCH_PATH"
    echo "    $boost_version at $BOOST_ARCH_PATH"
fi

if [ -d "$CGAL_ARCH_PATH" ]
then
    _foamAddLib $CGAL_ARCH_PATH/lib
fi

if [ -d "$BOOST_ARCH_PATH" ]
then
    _foamAddLib $BOOST_ARCH_PATH/lib
fi

unset boost_version cgal_version

# -----------------------------------------------------------------------------
