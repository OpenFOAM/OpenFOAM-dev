/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::Function1s::nonUniformTable

Description
    Non-uniform tabulated property function that linearly interpolates between
    the values.

    To speed-up the search of the non-uniform table a uniform jump-table is
    created on construction which is used for fast indirect addressing into
    the table.

Usage
    \nonUniformTable
        Property    | Description
        values      | List of value pairs
    \endnonUniformTable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    nonUniformTable;

        values
        (
            (280 999.87)
            (300 995.1)
            (350 973.7)
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef nonUniformTable1_H
#define nonUniformTable1_H

#include "Function1.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                       Class nonUniformTable Declaration
\*---------------------------------------------------------------------------*/

class nonUniformTable
:
    public FieldFunction1<scalar, nonUniformTable>
{
    // Private member data

        //- Lowest value in the nonUniformTable
        scalar low_;

        //- Highest value in the nonUniformTable
        scalar high_;

        //- Table values
        List<Tuple2<scalar, scalar>> values_;

        //- Increment derived from low_, high_ and values_.size()
        scalar delta_;

        List<label> jumpTable_;


protected:

        //- Return the lower index of the interval in the table
        //  corresponding to the given value
        inline label index(scalar x) const;


public:

    //- Runtime type information
    TypeName("nonUniformTable");


    // Constructors

        //- Construct from name and dictionary
        nonUniformTable(const word& name, const dictionary& dict);


    // Member Functions

        //- Return the non-uniform table of values
        const List<Tuple2<scalar, scalar>>& values() const
        {
            return values_;
        }

        //- Evaluate the function and return the result
        virtual scalar value(scalar x) const;

        //- Integrate between two scalar values
        virtual scalar integral(const scalar x1, const scalar x2) const;

        //- Evaluate the derivative of the function and return the result
        scalar dfdT(scalar T) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nonUniformTable1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
