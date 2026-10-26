/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::fallOffFunctions::Lindemann

Description
    Lindemann fall-off function

SourceFiles
    Lindemann_fallOffFunctionI.H

\*---------------------------------------------------------------------------*/

#ifndef Lindemann_fallOffFunction_H
#define Lindemann_fallOffFunction_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fallOffFunctions
{

/*---------------------------------------------------------------------------*\
                  Class Lindemann Declaration
\*---------------------------------------------------------------------------*/

class Lindemann
{

public:

    // Constructors

        //- Construct null
        inline Lindemann();

        //- Construct from dictionary
        inline Lindemann(const dictionary& dict);


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Lindemann";
        }

        //- Value
        inline scalar F(const scalar T, const scalar Pr) const;

        //- Derivative w.r.t. temperature
        inline scalar ddT
        (
            const scalar T,
            const scalar Pr,
            const scalar F
        ) const;

        //- Derivative w.r.t. reduced pressure
        inline scalar ddPr
        (
            const scalar T,
            const scalar Pr,
            const scalar F
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fallOffFunctions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
Ostream& operator<<(Ostream&, const fallOffFunctions::Lindemann&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Lindemann_fallOffFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
