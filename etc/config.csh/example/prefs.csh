#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     | Website:  https://openfoam.org
#   \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
#     config.csh/example/prefs.csh
#
# Description
#     Preset variables for the OpenFOAM configuration - C-Shell shell syntax.
#
#     The prefs.csh file will be sourced by the OpenFOAM etc/cshrc when it is
#     found by foamEtcFile.
#
# See also
#     'foamEtcFile -help' or 'foamEtcFile -list' for information about the
#     paths searched
#
#------------------------------------------------------------------------------

## Specify OpenFOAM ThirdParty compiler
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setenv WM_COMPILER_TYPE ThirdParty

## Specify compiler type
## ~~~~~~~~~~~~~~~~~~~~~
#setenv WM_COMPILER Clang

## Specify system openmpi
## ~~~~~~~~~~~~~~~~~~~~~~
#setenv WM_MPLIB SYSTEMOPENMPI

## Specify OpenFOAM ThirdParty openmpi with version
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setenv WM_MPLIB OPENMPI
#setenv OPENMPI_VERSION 2.1.1

## Specify options for decomposition libraries
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setenv SCOTCH_TYPE system
#setenv METIS_TYPE ThirdParty
#setenv METIS_VERSION 5.1.0
#setenv PARMETIS_TYPE ThirdParty
#setenv ZOLTAN_TYPE ThirdParty

## Specify system ParaView
#setenv ParaView_TYPE system

#------------------------------------------------------------------------------
