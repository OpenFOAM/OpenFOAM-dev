/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::Function1s::integratedNonUniformTable

Description
    Non-uniform tabulated property function that linearly interpolates between
    the values.

    To speed-up the search of the non-uniform table a uniform jump-table is
    created on construction which is used for fast indirect addressing into
    the table.

Usage
    \table
        Property    | Description
        values      | List of value pairs
    \endtable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    integratedNonUniformTable;

        values
        (
            (280 999.87)
            (300 995.1)
            (350 973.7)
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef integratedNonUniformTable1_H
#define integratedNonUniformTable1_H

#include "NonUniformTable1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                  Class integratedNonUniformTable Declaration
\*---------------------------------------------------------------------------*/

class integratedNonUniformTable
:
    public NonUniformTable<scalar>
{
    // Private member data

        List<scalar> intf_;

        List<scalar> intfByT_;


public:

    //- Runtime type information
    TypeName("integratedNonUniformTable");


    // Constructors

        //- Construct from entry name and dictionary
        integratedNonUniformTable(const word& name, const dictionary& dict);


    // Member Functions

        //- Integrate the function and return the result
        scalar intfdT(scalar T) const;

        //- Integrate the function divided by T and return the result
        scalar intfByTdT(scalar T) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
