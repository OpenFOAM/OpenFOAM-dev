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
    Foam::fallOffFunctions::Troe

Description
    The Troe fall-off function

SourceFiles
    Troe_fallOffFunctionI.H

\*---------------------------------------------------------------------------*/

#ifndef Troe_fallOffFunction_H
#define Troe_fallOffFunction_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fallOffFunctions
{

/*---------------------------------------------------------------------------*\
                     Class Troe Declaration
\*---------------------------------------------------------------------------*/

class Troe
{
    // Private Data

        scalar alpha_;

        scalar Tsss_, Ts_, Tss_;


public:

    // Constructors

        //- Construct from components
        inline Troe
        (
            const scalar alpha,
            const scalar Tsss,
            const scalar Ts,
            const scalar Tss
        );

        //- Construct from dictionary
        inline Troe(const dictionary& dict);


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Troe";
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
Ostream& operator<<(Ostream&, const fallOffFunctions::Troe&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Troe_fallOffFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
