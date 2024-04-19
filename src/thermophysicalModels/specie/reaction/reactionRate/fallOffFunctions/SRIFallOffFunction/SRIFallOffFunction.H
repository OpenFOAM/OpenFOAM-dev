/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::SRIFallOffFunction

Description
    The SRI fall-off function

SourceFiles
    SRIFallOffFunctionI.H

\*---------------------------------------------------------------------------*/

#ifndef SRIFallOffFunction_H
#define SRIFallOffFunction_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class SRIFallOffFunction;
Ostream& operator<<(Ostream&, const SRIFallOffFunction&);


/*---------------------------------------------------------------------------*\
                     Class SRIFallOffFunction Declaration
\*---------------------------------------------------------------------------*/

class SRIFallOffFunction
{
    // Private Data

        scalar a_, b_, c_, d_, e_;


public:

    // Constructors

        //- Construct from components
        inline SRIFallOffFunction
        (
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d,
            const scalar e
        );

        //- Construct from dictionary
        inline SRIFallOffFunction(const dictionary& dict);


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "SRI";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar Pr
        ) const;

        inline scalar ddT
        (
            const scalar T,
            const scalar Pr,
            const scalar F
        ) const;

        inline scalar ddPr
        (
            const scalar T,
            const scalar Pr,
            const scalar F
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const SRIFallOffFunction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SRIFallOffFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
