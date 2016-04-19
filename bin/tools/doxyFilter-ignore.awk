#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
# Script
#     doxyFilter-ignore.awk
#
# Description
#     - Prefix file contents with doxygen @file tag and %filePath% tag
#       that will be changed in a subsequent sed script
#     - Surround the contents of an entire file with @cond / @endcond
#       to skip documenting all classes/variables
#
#------------------------------------------------------------------------------
BEGIN {
   print "//! @file %filePath%"
   print "//! @cond OpenFOAMIgnoreAppDoxygen"
}

{ print }

END {
   print "//! @endcond"
}
#------------------------------------------------------------------------------
