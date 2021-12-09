/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "axisRotationMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(axisRotationMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        axisRotationMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::axisRotationMotion::axisRotationMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::axisRotationMotion::~axisRotationMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::axisRotationMotion::transformation() const
{
    const scalar t = time_.value() - time_.beginTime().value();

    if (t > small)
    {
        // Rotational position (in radians)
        const vector omega
        (
            t*degToRad(radialVelocity_.x()),
            t*degToRad(radialVelocity_.y()),
            t*degToRad(radialVelocity_.z())
        );

        const scalar magOmega = mag(omega);
        const quaternion R(omega/magOmega, magOmega);
        const septernion TR(septernion(-origin_)*R*septernion(origin_));

        DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

        return TR;
    }
    else
    {
        return septernion::I;
    }
}


bool Foam::solidBodyMotionFunctions::axisRotationMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("origin") >> origin_;
    SBMFCoeffs_.lookup("radialVelocity") >> radialVelocity_;

    return true;
}


// ************************************************************************* //
