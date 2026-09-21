/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "relativeMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(relativeMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        relativeMotion,
        PtrListDictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::relativeMotion::relativeMotion
(
    const word& name,
    const PtrListDictionary<solidBodyMotionFunction>& SBMFs,
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(name, SBMFCoeffs, runTime),
    SBMFs_(SBMFs),
    referenceName_(SBMFCoeffs.lookup<word>("relativeTo")),
    SBMF_(solidBodyMotionFunction::New(SBMFCoeffs, time_, typeName))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::relativeMotion::~relativeMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::relativeMotion::transformation() const
{
    // Get the reference transformation
    septernion TR = SBMFs_[referenceName_].transformation();

    // Set this body transformation relative to the reference
    TR *= SBMF_->transformation();

    DebugInFunction
        << "Time = " << time_.value() << " transformation: " << TR << endl;

    return TR;
}


// ************************************************************************* //
