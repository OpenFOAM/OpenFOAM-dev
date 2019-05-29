/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::Function1Types::quadraticRamp

Description
    Quadratic ramp function starting from 0 and increasing quadraticRampally
    to 1 from \c t_0 over the \c duration and remaining at 1 thereafter.

See also
    Foam::Function1Types::ramp

SourceFiles
    quadraticRamp.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticRamp_H
#define quadraticRamp_H

#include "ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class quadraticRamp Declaration
\*---------------------------------------------------------------------------*/

class quadraticRamp
:
    public ramp
{

public:

    // Runtime type information
    TypeName("quadraticRamp");


    // Constructors

        //- Construct from entry name and dictionary
        quadraticRamp
        (
            const word& entryName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~quadraticRamp();


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const quadraticRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quadraticRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
