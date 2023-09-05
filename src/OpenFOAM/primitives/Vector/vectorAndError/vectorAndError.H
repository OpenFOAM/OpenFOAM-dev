/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::vectorAndError

Description
    Vector of scalarAndError-s.

\*---------------------------------------------------------------------------*/

#ifndef vectorAndError_H
#define vectorAndError_H

#include "scalarAndError.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class vectorAndError Declaration
\*---------------------------------------------------------------------------*/

class vectorAndError
:
    public Vector<scalarAndError>
{
public:

    // Constructors

        //- Inherit Constructors
        using Vector<scalarAndError>::Vector;


    // Member Functions

        //- Return the value
        inline Vector<scalar> value() const
        {
            return Vector<scalar>(x().value, y().value, z().value);
        }

        //- Return the error
        inline Vector<scalar> error() const
        {
            return Vector<scalar>(x().error, y().error, z().error);
        }
};


//- Data associated with the vector-and-error type is contiguous
template<>
inline bool contiguous<vectorAndError>() { return true; }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
