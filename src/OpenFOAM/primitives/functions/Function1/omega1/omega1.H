/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::Function1s::omega

Description
    Convenience class to handle the input of time-varying rotational speed.
    Reads an \c omega Function1 entry with default units of [rad/s]. For
    backwards compatibility this will also alternatively read an \c rpm entry
    with default units of [rpm].

Usage
    For specifying a constant rotational speed of an MRF zone:
    \verbatim
        MRF
        {
            cellZone    rotor;

            origin     (0 0 0);
            axis       (0 0 1);

            omega      6.28319;  // <-- Basic specification in [rad/s]

        //  omega      60 [rpm]; // <-- Equivalent specification with unit
                                 //     conversion from [rpm]

        //  rpm        60;       // <-- Equivalent backwards compatible
                                 //     specification for rpm. May be
                                 //     removed in future.
        }
    \endverbatim
    or for a tabulated ramped rotational speed of a solid body:
    \verbatim
        mover
        {
            type            motionSolver;

            libs            ("libfvMotionSolvers.so");

            motionSolver    solidBody;

            cellZone        innerCylinder;

            solidBodyMotionFunction  rotatingMotion;

            origin      (0 0 0);
            axis        (0 1 0);

            omega       table            // <-- Basic specification in [rad/s]
            (
                (0      0)
                (0.01   628.319)
                (0.022  628.319)
                (0.03   418.879)
                (100    418.879)
            );

        //  omega                        // <-- Equivalent specification with
        //  {                            //     unit conversion from [rpm].
        //      type        table;       //     Note that the dictionary form
        //      units       ([s] [rpm]); //     is needed for the units entry.
        //      values
        //      (
        //          (0      0)
        //          (0.01   6000)
        //          (0.022  6000)
        //          (0.03   4000)
        //          (100    4000)
        //      );
        //  }

        //  rpm         table            // <-- Equivalent backwards compatible
        //  (                            //     specification for rpm. May be
        //      (0      0)               //     removed in future.
        //      (0.01   6000)
        //      (0.022  6000)
        //      (0.03   4000)
        //      (100    4000)
        //  );
        }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    omega1.C
    omega1I.H

\*---------------------------------------------------------------------------*/

#ifndef omega1_H
#define omega1_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class omega Declaration
\*---------------------------------------------------------------------------*/

class omega
{
    // Private Data

        //- Reference to the time units
        const unitConversion& tUnits_;

        //- The omega function
        const autoPtr<Function1<scalar>> omega_;


    // Private Member Functions

        //- Read and construct the omega function
        autoPtr<Function1<scalar>> init
        (
            const Time& time,
            const dictionary& dict
        );


public:

    // Constructors

        //- Construct from dictionary
        omega(const Time& time, const dictionary& dict);

        //- Copy constructor
        omega(const omega&);

        //- Construct and return a clone
        autoPtr<omega> clone() const
        {
            return autoPtr<omega>(new omega(*this));
        }


    // Member Functions

        //- Return value for time t
        inline scalar value(const scalar t) const;

        //- Return the integral between times t1 and t2
        inline scalar integral
        (
            const scalar t1,
            const scalar t2
        ) const;

        //- Write data to dictionary stream
        void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const omega&) = delete;
};


void writeEntry(Ostream& os, const omega& a);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "omega1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
