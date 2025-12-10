/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "multiValveEnginePistonMotion.H"
#include "multiValveEngine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(multiValveEnginePistonMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        multiValveEnginePistonMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiValveEnginePistonMotion::
multiValveEnginePistonMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    fluidRegionName_()
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiValveEnginePistonMotion::
~multiValveEnginePistonMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::multiValveEnginePistonMotion::
transformation() const
{
    const fvMesh& fluidMesh = time_.lookupObject<fvMesh>(fluidRegionName_);

    const fvMeshMover& fluidMover = fluidMesh.mover();

    if (!isA<fvMeshMovers::multiValveEngine>(fluidMover))
    {
        FatalErrorInFunction
            << "The " << typeName << " motion function requires the motion"
            << " solver for the associated fluid region (i.e., "
            << fluidRegionName_ << ") to be of type "
            << fvMeshMovers::multiValveEngine::typeName
            << exit(FatalError);
    }

    const fvMeshMovers::multiValveEngine& fluidEngineMover =
        refCast<const fvMeshMovers::multiValveEngine>(fluidMover);

    return
        septernion
        (
            fluidEngineMover.piston.position()
           *fluidEngineMover.piston.axis
        );
}


bool Foam::solidBodyMotionFunctions::multiValveEnginePistonMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    fluidRegionName_ = SBMFCoeffs.lookup<word>("fluidRegion");

    return true;
}


// ************************************************************************* //
