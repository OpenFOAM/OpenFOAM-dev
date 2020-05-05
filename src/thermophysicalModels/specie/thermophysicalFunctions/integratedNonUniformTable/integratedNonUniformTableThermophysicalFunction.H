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
    Foam::thermophysicalFunctions::integratedNonUniformTable

Description
    Non-uniform tabulated property function that linearly interpolates between
    the values.

    To speed-up the search of the non-uniform table a uniform jump-table is
    created on construction which is used for fast indirect addressing into
    the table.

Usage
    \integratedNonUniformTable
        Property    | Description
        values      | List of (temperature property) value pairs
    \endintegratedNonUniformTable

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

#ifndef integratedNonUniformTableThermophysicalFunction_H
#define integratedNonUniformTableThermophysicalFunction_H

#include "nonUniformTableThermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{

/*---------------------------------------------------------------------------*\
                  Class integratedNonUniformTable Declaration
\*---------------------------------------------------------------------------*/

class integratedNonUniformTable
:
    public nonUniformTable
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

        //- Construct from dictionary
        integratedNonUniformTable(const dictionary& dict);


    // Member Functions

        //- Integrate the function and return the result
        scalar intfdT(scalar p, scalar T) const;

        //- Integrate the function divided by T and return the result
        scalar intfByTdT(scalar p, scalar T) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermophysicalFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
