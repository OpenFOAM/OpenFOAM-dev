/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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
    Foam::Function1s::AntoineExtended

Description
    Extended Antoine equation, used to describe vapour pressure as a function
    of temperature

    \f[
        \log (p) = A + \frac{B}{C + T} + D T + E \log (T) + F T^G
    \f]

    Note that unlike the non-extended Antoine equation, this uses natural
    logarithms and the second term has a positive sign

SourceFiles
    AntoineExtended.C

\*---------------------------------------------------------------------------*/

#ifndef AntoineExtended_H
#define AntoineExtended_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class AntoineExtended Declaration
\*---------------------------------------------------------------------------*/

class AntoineExtended
:
    public FieldFunction1<scalar, AntoineExtended>
{
    // Private data

        //- Coefficients
        const scalar A_, B_, C_, D_, E_, G_, F_;


public:

    //- Runtime type information
    TypeName("AntoineExtended");


    // Constructors

        //- Construct from a name and dictionary
        AntoineExtended
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );


    //- Destructor
    virtual ~AntoineExtended();


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual inline scalar value(const scalar T) const;

        //- Return the derivative at a scalar coordinate
        virtual inline scalar derivative(const scalar T) const;

        //- Return the integral between two scalar coordinates
        virtual inline scalar integral(const scalar T1, const scalar T2) const;

        //- Write in dictionary format
        virtual void write(Ostream& os, const unitSets& units) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AntoineExtendedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
