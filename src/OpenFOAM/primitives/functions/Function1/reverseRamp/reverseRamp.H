/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::Function1s::reverseRamp

Description
    A wrapper class to reverse any ramp function such that the result starts
    from 1 decreasing to 0 from \c start over the \c duration and remaining at 0
    thereafter.

    Usage for scaling a vector:
    \verbatim
        <name>
        {
            scale
            {
                type        reverseRamp;

                ramp        linearRamp;

                start       0;
                duration    10;
            }

            value
            {
                type        sine;

                frequency   10;
                amplitude   1;
                scale       (1 0.1 0);
                level       (10 1 0);
            }
        }
    \endverbatim

See also
    Foam::Function1s::Ramp

SourceFiles
    reverseRamp.C

\*---------------------------------------------------------------------------*/

#ifndef reverseRamp_H
#define reverseRamp_H

#include "Ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class reverseRamp Declaration
\*---------------------------------------------------------------------------*/

class reverseRamp
:
    public Ramp<reverseRamp>
{
    // Private Data

        //- Standard ramp function to reverse
        autoPtr<Function1<scalar>> ramp_;


public:

    // Runtime type information
    TypeName("reverseRamp");


    // Constructors

        //- Construct from name and dictionary
        reverseRamp
        (
            const word& name,
            const dictionary& dict
        );

        //- Copy constructor
        reverseRamp(const reverseRamp& rr);


    //- Destructor
    virtual ~reverseRamp();


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;

        //- Return the integral between times t1 and t2
        virtual inline scalar integral
        (
            const scalar t1,
            const scalar t2
        ) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reverseRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reverseRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
