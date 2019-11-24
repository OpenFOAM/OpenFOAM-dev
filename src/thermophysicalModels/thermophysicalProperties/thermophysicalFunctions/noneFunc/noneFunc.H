/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::noneFunc

Description
    Undefined function which returns an error when called.

\*---------------------------------------------------------------------------*/

#ifndef noneFunc_H
#define noneFunc_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class noneFunc;
Ostream& operator<<(Ostream& os, const noneFunc& f);


/*---------------------------------------------------------------------------*\
                           Class noneFunc Declaration
\*---------------------------------------------------------------------------*/

class noneFunc
:
    public thermophysicalFunction
{
    // Private member data

        //- Name of dictionary from which this function is instantiated
        fileName dictName_;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        noneFunc(const dictionary& dict);


    // Member Functions

        //- Evaluate the function and return the result
        scalar f(scalar p, scalar T) const;

        //- Write the function coefficients
        void writeData(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const noneFunc& f);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
