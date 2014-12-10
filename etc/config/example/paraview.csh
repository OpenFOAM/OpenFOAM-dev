#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#     config/example/paraview.csh
#
# Description
#     Example of chaining to the standard config/paraview.csh with a
#     different ParaView_VERSION
#
# Note
#     This file could be copied to a user or site location, but should never
#     replace the default shipped version as this will cause an infinite loop
#
#------------------------------------------------------------------------------

#
# Use other (shipped) paraview.csh with a different ParaView_VERSION
#

set foamFile=`$WM_PROJECT_DIR/bin/foamEtcFile -mode o config/paraview.csh`
if ( $status == 0 ) source $foamFile ParaView_VERSION=3.12.0

unset foamFile

# -----------------------------------------------------------------------------
