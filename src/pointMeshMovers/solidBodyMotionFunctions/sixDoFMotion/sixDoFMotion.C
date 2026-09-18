/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "sixDoFMotion.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(sixDoFMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        sixDoFMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::sixDoFMotion::sixDoFMotion
(
    const word& name,
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(name, SBMFCoeffs, runTime),
    CofG_(SBMFCoeffs_.lookup("CofG")),
    translation_
    (
        Function1<vector>::New
        (
            "translation",
            time_.userUnits(),
            dimensions::length,
            SBMFCoeffs_
        ).ptr()
    ),
    rotation_
    (
        Function1<vector>::New
        (
            "rotation",
            time_.userUnits(),
            units::degrees,
            SBMFCoeffs_
        ).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::sixDoFMotion::~sixDoFMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::sixDoFMotion::transformation() const
{
    const scalar t = time_.value();

    const vector translation = translation_->value(t);
    const vector rotation = rotation_->value(t);

    const quaternion R(quaternion::XYZ, rotation);
    const septernion TR(septernion(-CofG_ - translation)*R*septernion(CofG_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


// ************************************************************************* //
