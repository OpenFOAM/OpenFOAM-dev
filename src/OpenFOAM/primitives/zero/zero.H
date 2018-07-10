/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::zero

Description
    A class representing the concept of 0 used to avoid unnecessary
    manipulations for objects that are known to be zero at compile-time.

SourceFiles
    zeroI.H

\*---------------------------------------------------------------------------*/

#ifndef zero_H
#define zero_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class zero Declaration
\*---------------------------------------------------------------------------*/

class zero
{
public:

    typedef zero value_type;

    // Constructors

        //- Construct null
        zero()
        {}


    // Member operators

        //- Return 0 for bool
        inline operator bool() const
        {
            return 0;
        }

        //- Return 0 for label
        inline operator label() const
        {
            return 0;
        }

        //- Return 0 for float
        inline operator float() const
        {
            return 0;
        }

        //- Return 0 for double
        inline operator double() const
        {
            return 0;
        }

        //- Return 0 for double
        inline operator long double() const
        {
            return 0;
        }
};


// Global zero
static const zero Zero;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "zeroI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
