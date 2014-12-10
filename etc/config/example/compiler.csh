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
#     config/example/compiler.csh
#
# Description
#     Example of fine tuning ThirdParty compiler settings for OpenFOAM
#     Sourced from OpenFOAM-<VERSION>/etc/config/settings.csh
#
#------------------------------------------------------------------------------

# Modified compiler settings
switch ("$WM_COMPILER")
case Gcc46:
case Gcc46++0x:
    set gcc_version=gcc-4.6.0
    set gmp_version=gmp-5.0.1
    set mpfr_version=mpfr-2.4.2
    set mpc_version=mpc-0.8.1
    breaksw
case Gcc45:
case Gcc45++0x:
    set gcc_version=gcc-4.5.2
    set gmp_version=gmp-5.0.1
    set mpfr_version=mpfr-2.4.2
    set mpc_version=mpc-0.8.1
    breaksw
endsw

# ----------------------------------------------------------------- end-of-file
