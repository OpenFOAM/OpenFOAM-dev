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
    Foam::scalarAndError

Description
    Class to encapsulate a scalar value and an associated round-off error.
    The error is tracked through operations performed between scalarAndError
    variables.

SourceFiles
    scalarAndErrorI.H

\*---------------------------------------------------------------------------*/

#ifndef scalarAndError_H
#define scalarAndError_H

#include "scalar.H"
#include "zero.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scalarAndError Declaration
\*---------------------------------------------------------------------------*/

class scalarAndError
{
public:

    //- Component type
    typedef scalarAndError cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of scalarAndError is 0
        static const direction rank = 0;

        //- Number of components in scalarAndError is 1
        static const direction nComponents = 1;


    // Public Data

        //- The value
        scalar value;

        //- The error
        scalar error;


    // Constructors

        //- Default construct
        inline scalarAndError()
        {}

        //- Construct zero
        inline scalarAndError(const zero&)
        :
            value(scalar(0)),
            error(scalar(0))
        {}

        //- Construct from a scalar
        inline scalarAndError(const scalar v)
        :
            value(v),
            error(small*mag(v))
        {}

        //- Construct from a scalar and an error
        inline scalarAndError(const scalar v, const scalar e)
        :
            value(v),
            error(e)
        {}


    // Member Operators

        //- Negate
        inline scalarAndError operator-() const
        {
            return scalarAndError(- value, error);
        }

        //- Add another scalar-and-error to this one
        inline void operator+=(const scalarAndError& sae)
        {
            value += sae.value;
            error += mag(sae.error);
        }

        //- Subtract another scalar-and-error from this one
        inline void operator-=(const scalarAndError& sae)
        {
            this->operator+=(- sae);
        }

        //- Multiply another scalar-and-error with this one
        inline void operator*=(const scalarAndError& sae)
        {
            error = mag(value*error) + mag(sae.value*sae.error);
            value *= sae.value;
        }

        //- Divide this scalar-and-error by another
        inline void operator/=(const scalarAndError& sae)
        {
            error = mag(error/sae.value) + mag(sae.error/sqr(sae.value));
            value /= sae.value;
        }
};


//- Data associated with the scalar-and-error type is contiguous
template<>
inline bool contiguous<scalarAndError>() { return true; }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarAndErrorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
