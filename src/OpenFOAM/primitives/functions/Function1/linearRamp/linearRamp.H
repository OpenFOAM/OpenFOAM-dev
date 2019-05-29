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
    Foam::Function1Types::linearRamp

Description
    Linear ramp function starting from 0 and increasing linearRamply to 1 from
    \c start over the \c duration and remaining at 1 thereafter.

See also
    Foam::Function1Types::ramp

SourceFiles
    linearRamp.C

\*---------------------------------------------------------------------------*/

#ifndef linearRamp_H
#define linearRamp_H

#include "ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class linearRamp Declaration
\*---------------------------------------------------------------------------*/

class linearRamp
:
    public ramp
{

public:

    // Runtime type information
    TypeName("linearRamp");


    // Constructors

        //- Construct from entry name and dictionary
        linearRamp
        (
            const word& entryName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearRamp();


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "linearRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
