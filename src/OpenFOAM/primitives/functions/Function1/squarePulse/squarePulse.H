/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::Function1s::squarePulse

Description
    Square pulse function. Takes a value of 0 before \c start, then rises
    instantaneously to 1, remains at 1 for \c duration, then drops
    instantaneously back to 0.

    Usage:
    \verbatim
        <name> squarePulse;
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
            type      squarePulse;
            start     10;
            duration  20;
        }
    \endverbatim

    Where:
    \table
        Property | Description  | Required | Default value
        start    | Start time   | no       | 0
        duration | Duration     | yes      |
    \endtable

See also
    Foam::Function1

SourceFiles
    squarePulse.C

\*---------------------------------------------------------------------------*/

#ifndef squarePulse_H
#define squarePulse_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class squarePulse Declaration
\*---------------------------------------------------------------------------*/

class squarePulse
:
    public FieldFunction1<scalar, squarePulse>
{
    // Private data

        //- Start-time of the ramp function
        scalar start_;

        //- Duration of the ramp function
        scalar duration_;


private:

    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& dict);


public:

    // Runtime type information
    TypeName("squarePulse");


    // Constructors

        //- Construct from name and dictionary
        squarePulse
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~squarePulse();


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
        void operator=(const squarePulse&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "squarePulseI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
