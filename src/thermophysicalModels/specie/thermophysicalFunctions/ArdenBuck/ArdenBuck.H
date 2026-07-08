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
    Foam::Function1s::ArdenBuck

Description
    Arden-Buck equation for the vapour pressure of moist air

SourceFiles
    ArdenBuck.C

\*---------------------------------------------------------------------------*/

#ifndef ArdenBuck_H
#define ArdenBuck_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                           Class ArdenBuck Declaration
\*---------------------------------------------------------------------------*/

class ArdenBuck
:
    public FieldFunction1<scalar, ArdenBuck>
{
    // Private Static Member Data

        //- Zero degrees Celcius, in whatever the standard temperature units are
        static const scalar zeroC_;

        //- Coefficients
        static const scalar A_, B_, C_, D_;


    // Private Member Functions

        //- Exponent divided by the temperature
        inline scalar xByTC(const scalar TC) const;


public:

    //- Runtime type information
    TypeName("ArdenBuck");


    // Constructors

        //- Construct from a name and dictionary
        ArdenBuck
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );

        //- Construct from name and Istream
        ArdenBuck(const word& name, const unitSets& units, Istream& is);


    //- Destructor
    virtual ~ArdenBuck();


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

#include "ArdenBuckI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
