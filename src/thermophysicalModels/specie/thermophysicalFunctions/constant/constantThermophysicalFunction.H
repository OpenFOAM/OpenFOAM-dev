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
    Foam::thermophysicalFunctions::constant

Description
    Constant property function.

Usage
    \constant
        Property    | Description
        value       | Constant property value
    \endconstant

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    constant;
        value   999.87;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef constantThermophysicalFunction_H
#define constantThermophysicalFunction_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{

/*---------------------------------------------------------------------------*\
                          Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public thermophysicalFunction
{
    // Private member data

        //- Constant value
        scalar value_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary
        constant(const dictionary& dict);


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
