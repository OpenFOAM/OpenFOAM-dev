/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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
    Foam::Function1s::laminarBL

Description
    Laminar boundary-layer function yHat*(2 - yHat)

    where yHat is the normalised distance to the wall.

    Usage:
    \verbatim
        <name> laminarBL;
    \endverbatim

See also
    Foam::Function1
    Foam::Function1s::turbulentBL

SourceFiles
    laminarBL.C

\*---------------------------------------------------------------------------*/

#ifndef laminarBL_H
#define laminarBL_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class laminarBL Declaration
\*---------------------------------------------------------------------------*/

class laminarBL
:
    public FieldFunction1<scalar, laminarBL>
{

public:

    // Runtime type information
    TypeName("laminarBL");


    // Constructors

        //- Construct from name and dictionary
        laminarBL
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );


    //- Destructor
    virtual ~laminarBL();


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual inline scalar value(const scalar x) const;

        //- Return the derivative at a scalar coordinate
        virtual inline scalar derivative(const scalar x) const;

        //- Return the integral between two scalar coordinates
        virtual inline scalar integral(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void write(Ostream& os, const unitSets& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laminarBL&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "laminarBLI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
