/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::one

Description
    A class representing the concept of 1 (scalar(1)) used to avoid
    unnecessary manipulations for objects which are known to be one at
    compile-time.

SourceFiles
    oneI.H

\*---------------------------------------------------------------------------*/

#ifndef one_H
#define one_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class one Declaration
\*---------------------------------------------------------------------------*/

class one
{
public:

    typedef one value_type;

    // Constructors

        //- Construct null
        one()
        {}


    // Member Operators

        //- Return 1 for label
        inline operator label() const
        {
            return 1;
        }

        //- Return 1 for float
        inline operator float() const
        {
            return 1;
        }

        //- Return 1 for double
        inline operator double() const
        {
            return 1;
        }

        //- Return 1 for double
        inline operator long double() const
        {
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "oneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
