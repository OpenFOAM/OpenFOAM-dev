#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     | Website:  https://openfoam.org
#   \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#     doxyFilter.awk
#
# Description
#     Converts cocoon style sentinel strings into doxygen style strings
#
#     Assumes comment strings are formatted as follows
#         //- General description
#         //  more information
#         //  and even more information
#     This should be re-formatted as the following
#         //! general description
#         /*!
#         more information
#         and even more information
#         */
#     The intermediate "/*! ... */" block is left-justified to handle
#     possible verbatim text
#------------------------------------------------------------------------------

BEGIN {
    state = 0
}

/^ *\/\/-/ {
    state = 1
    sub(/\/\/-/, "//!")
    print
    next
}

/^ *\/\// {
    # Start comment block
    if (state == 1)
    {
        printf "/*! "
        state = 2
    }

    # Inside comment block
    if (state == 2)
    {
        if (!sub(/^ *\/\/  /, ""))
        {
            sub(/^ *\/\//, "")
        }
    }

    print
    next
}

{
    # End comment block
    if (state == 2)
    {
        printf "*/ "
    }
    state = 0
    print
    next
}

#------------------------------------------------------------------------------
