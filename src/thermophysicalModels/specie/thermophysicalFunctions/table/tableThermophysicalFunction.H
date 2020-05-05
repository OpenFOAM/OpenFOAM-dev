/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::thermophysicalFunctions::table

Description
    Tabulated property function that linearly interpolates between
    the table values.

Usage
    \table
        Property    | Description
        Tlow        | Lower temperature limit of the table
        Thigh       | Upper temperature limit of the table
        values      | Property values, assumed uniformly distributed
    \endtable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    table;

        Tlow    280;
        Thigh   350;

        values
        (
            999.87
            991.53
            973.7
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef tableThermophysicalFunction_H
#define tableThermophysicalFunction_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{

/*---------------------------------------------------------------------------*\
                           Class table Declaration
\*---------------------------------------------------------------------------*/

class table
:
    public thermophysicalFunction
{
    // Private member data

        //- Name of dictionary from which this function is instantiated
        fileName dictName_;

        //- Lowest temperature in the table
        scalar Tlow_;

        //- Highest temperature in the table
        scalar Thigh_;

        //- Table values
        List<scalar> values_;

        //- Temperature increment derived from Tlow_, Thigh_ and values_.size()
        scalar deltaT_;


public:

    //- Runtime type information
    TypeName("table");


    // Constructors

        //- Construct from dictionary
        table(const dictionary& dict);


    // Member Functions

        //- Evaluate the function and return the result
        scalar f(scalar p, scalar T) const;

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermophysicalFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
