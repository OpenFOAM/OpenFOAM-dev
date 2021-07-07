/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
    Foam::Function1s::UniformTable

Description
    Tabulated property function that linearly interpolates between
    the UniformTable values.

Usage
    \table
        Property    | Description
        low         | Lower limit of the table
        high        | Upper limit of the table
        values      | Values, assumed uniformly distributed
    \endtable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    uniformTable;

        low    280;
        high   350;

        values
        (
            999.87
            991.53
            973.7
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef UniformTable1_H
#define UniformTable1_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class UniformTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformTable
:
    public FieldFunction1<Type, UniformTable<Type>>
{
    // Private member data

        //- Name of dictionary from which this function is instantiated
        fileName dictName_;

        //- Lowest value in the table
        scalar low_;

        //- Highest value in the table
        scalar high_;

        //- Table values
        List<Type> values_;

        //- Increment derived from low_, high_ and values_.size()
        scalar delta_;


public:

    //- Runtime type information
    TypeName("uniformTable");


    // Constructors

        //- Construct from name and dictionary
        UniformTable(const word& name, const dictionary& dict);


    // Member Functions

        //- Evaluate the function and return the result
        virtual Type value(scalar x) const;

        //- Integrate between two scalar values
        virtual Type integral(const scalar x1, const scalar x2) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UniformTable1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
