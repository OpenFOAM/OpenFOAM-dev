/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "mathematicalConstants.H"
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "None.H"
#include "Constant.H"
#include "Uniform.H"
#include "ZeroConstant.H"
#include "OneConstant.H"
#include "Polynomial1.H"
#include "Sine.H"
#include "Square.H"
#include "Table.H"
#include "UniformTable1.H"
#include "NonUniformTable1.H"
#include "EmbeddedTableReader.H"
#include "FoamTableReader.H"
#include "Scale.H"
#include "CodedFunction1.H"

typedef Foam::solidBodyMotionFunctions::sixDoFMotion::translationRotationVectors
trvType;

template<>
const char* const trvType::vsType::typeName = "vector2Vector";

template<>
const char* const trvType::vsType::componentNames[] = {"x", "y"};

template<>
const trvType trvType::vsType::vsType::zero
(
    trvType::uniform(vector::uniform(0))
);

template<>
const trvType trvType::vsType::one(trvType::uniform(vector::uniform(1)));

template<>
const trvType trvType::vsType::max(trvType::uniform(vector::uniform(vGreat)));

template<>
const trvType trvType::vsType::min(trvType::uniform(vector::uniform(-vGreat)));

template<>
const trvType trvType::vsType::rootMax
(
    trvType::uniform(vector::uniform(rootVGreat))
);

template<>
const trvType trvType::vsType::rootMin
(
    trvType::uniform(vector::uniform(-rootVGreat))
);

namespace Foam
{
    makeFunction1s(trvType);

    defineTableReader(trvType);
    makeTableReader(Embedded, trvType);
    makeTableReader(Foam, trvType);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::sixDoFMotion::sixDoFMotion
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

Foam::solidBodyMotionFunctions::sixDoFMotion::~sixDoFMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::sixDoFMotion::transformation() const
{
    const scalar t = time_.value();

    translationRotationVectors TRV = translationRotation_->value(t);

    // Convert the rotational motion from deg to rad
    TRV[1] *= pi/180.0;

    quaternion R(quaternion::XYZ, TRV[1]);
    septernion TR(septernion(-CofG_ + -TRV[0])*R*septernion(CofG_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::sixDoFMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    translationRotation_ = Function1<translationRotationVectors>::New
    (
        "translationRotation",
        SBMFCoeffs
    );

    SBMFCoeffs_.lookup("CofG") >> CofG_;

    return true;
}


// ************************************************************************* //
