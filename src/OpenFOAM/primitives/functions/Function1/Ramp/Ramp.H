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
    Foam::Function1s::Ramp

Description
    Ramp function base class for the set of scalar functions starting from 0 and
    increasing monotonically to 1 from \c start over the \c duration and
    remaining at 1 thereafter.

Usage
    \verbatim
        <name> <rampFunction>;
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
            type      <rampFunction>;
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
                type        <rampFunction>;

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
    Foam::Function1
    Foam::Function1s::Scale

SourceFiles
    Ramp.C

\*---------------------------------------------------------------------------*/

#ifndef Ramp_H
#define Ramp_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Ramp Declaration
\*---------------------------------------------------------------------------*/

template<class Function1Type>
class Ramp
:
    public FieldFunction1<scalar, Function1Type>
{
protected:

    // Protected data

        //- Start-time of the ramp function
        scalar start_;

        //- Duration of the ramp function
        scalar duration_;

        //- Return the end-time of the ramp function
        inline scalar end() const
        {
            return start_ + duration_;
        }

        //- Simple linear ramp function which form the basis of many more
        //  complex ramp functions
        inline scalar linearRamp(const scalar t) const
        {
            return max(min((t - start_)/duration_, 1), 0);
        }

        //- The gradient of the linear ramp within the bounds of the ramp
        //  duration
        inline scalar dLinearRampDt() const
        {
            return 1/duration_;
        }


private:

    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Constructors

        //- Construct from name and dictionary
        Ramp
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Ramp();


    // Member Functions

        //- Return value for time t
        virtual scalar value(const scalar t) const = 0;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Ramp<Function1Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Ramp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
