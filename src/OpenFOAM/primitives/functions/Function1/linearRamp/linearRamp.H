/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
    Foam::Function1s::linearRamp

Description
    Linear ramp function starting from 0 and increasing linearRamply to 1 from
    \c start over the \c duration and remaining at 1 thereafter.

Usage
    \verbatim
        <name> linearRamp;
        <name>Coeffs
        {
            start     10;
            duration  20;
        }
    \endverbatim
    or
    \verbatim
        <name>
        {
            type      linearRamp;
            start     10;
            duration  20;
        }
    \endverbatim
    or including a constant scaling factor:
    \verbatim
        <name>
        {
            type      scale;

            scale
            {
                type        linearRamp;

                start       0;
                duration    10;
            }

            value       100;
        }
    \endverbatim
    for more details of scaling the ramp see Foam::Function1s::Scale.

    Where:
    \table
        Property | Description  | Required | Default value
        start    | Start time   | no       | 0
        duration | Duration     | yes      |
    \endtable

See also
    Foam::Function1s::Ramp
    Foam::Function1s::Scale

SourceFiles
    linearRamp.C

\*---------------------------------------------------------------------------*/

#ifndef linearRamp_H
#define linearRamp_H

#include "Ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class linearRamp Declaration
\*---------------------------------------------------------------------------*/

class linearRamp
:
    public Ramp<linearRamp>
{

public:

    // Runtime type information
    TypeName("linearRamp");


    // Constructors

        //- Construct from name and dictionary
        linearRamp
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearRamp();


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;

        //- Return the integral between times t1 and t2
        virtual inline scalar integral
        (
            const scalar t1,
            const scalar t2
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "linearRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
