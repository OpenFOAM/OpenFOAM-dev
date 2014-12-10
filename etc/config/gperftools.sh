#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
#     config/gperftools.sh
#
# Description
#     Setup file for gperftools binaries libraries.
#
#------------------------------------------------------------------------------

version=svn
gperftools_install=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER

GPERFTOOLS_VERSION=gperftools-$version
GPERFTOOLS_ARCH_PATH=$gperftools_install/$GPERFTOOLS_VERSION

export PATH=$GPERFTOOLS_ARCH_PATH/bin:$PATH
export LD_LIBRARY_PATH=$GPERFTOOLS_ARCH_PATH/lib:$LD_LIBRARY_PATH

# -----------------------------------------------------------------------------
