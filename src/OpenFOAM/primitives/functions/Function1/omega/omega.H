/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    User convenience class to handle the input of time-varying rotational speed
    in rad/s if \c omega is specified or rpm if \c rpm is specified.

Usage
    For specifying the rotational speed in rpm of an MRF zone:
    \verbatim
        MRF
        {
            cellZone    rotor;

            origin     (0 0 0);
            axis       (0 0 1);

            rpm        60;
        }
    \endverbatim
    or the equivalent specified in rad/s:
    \verbatim
        MRF
        {
            cellZone    rotor;

            origin     (0 0 0);
            axis       (0 0 1);

            omega      6.28319;
        }
    \endverbatim
    or for a tabulated ramped rotational speed of a solid body:
    \verbatim
        mover
        {
            type            motionSolver;

            libs            ("libfvMeshMovers.so" "libfvMotionSolvers.so");

            motionSolver    solidBody;

            cellZone        innerCylinder;

            solidBodyMotionFunction  rotatingMotion;

            origin      (0 0 0);
            axis        (0 1 0);

            rpm         table
            (
                (0    0)
                (0.01  6000)
                (0.022  6000)
                (0.03  4000)
                (100   4000)
            );
        }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    omega.C

\*---------------------------------------------------------------------------*/

#ifndef omega_H
#define omega_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class omega Declaration
\*---------------------------------------------------------------------------*/

class omega
{
    // Private Data

        //- True if the input specification is rpm rather than omega
        bool rpm_;

        //- 1 for omega, pi/30 for rpm
        scalar omegaFactor_;

        //- The omega function
        autoPtr<Function1<scalar>> omega_;


public:

    // Constructors

        //- Construct null
        omega();

        //- Construct from dictionary
        omega(const dictionary& dict);

        //- Copy constructor
        omega(const omega&);


    // Member Functions

        //- Update omega function from given dictionary
        bool read(const dictionary& dict);

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

#include "omegaI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
