#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#     config/ensight.sh
#
# Description
#     Setup file for Ensight
#     Sourced from OpenFOAM-*/etc/bashrc
#
#------------------------------------------------------------------------------

# fallback value
if [ ! -d "$CEI_HOME" ]
then
    export CEI_HOME=/usr/local/ensight/CEI
fi

if [ -r $CEI_HOME ]
then

    # special treatment for 32bit OpenFOAM and 64bit Ensight
    if [ "$WM_ARCH" = linux -a `uname -m` = x86_64 ]
    then
        export CEI_ARCH=linux_2.6_32
    fi

    # add to path if required
    if [ "$CEI_HOME/bin/ensight" != "`which ensight 2>/dev/null`" ]
    then
        export PATH=$CEI_HOME/bin:$PATH
    fi

    export ENSIGHT9_INPUT=dummy
    export ENSIGHT9_READER=$FOAM_LIBBIN
else
    unset CEI_HOME
fi

# -----------------------------------------------------------------------------
