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

#include "multiValveEngineValveMotion.H"
#include "multiValveEngine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(multiValveEngineValveMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        multiValveEngineValveMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiValveEngineValveMotion::
multiValveEngineValveMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    fluidRegionName_(),
    valveName_(),
    valveIndex_(-1)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::multiValveEngineValveMotion::
~multiValveEngineValveMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::multiValveEngineValveMotion::
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

    if (valveIndex_ == -1)
    {
        forAll(fluidEngineMover.valves, valvei)
        {
            if (fluidEngineMover.valves[valvei].name == valveName_)
            {
                valveIndex_ = valvei;
                break;
            }
        }
    }

    if (valveIndex_ == -1)
    {
        wordList valveNames(fluidEngineMover.valves.size());
        forAll(fluidEngineMover.valves, valvei)
        {
            valveNames[valvei] = fluidEngineMover.valves[valvei].name;
        }

        FatalErrorInFunction
            << "A valve named '" << valveName_ << "' was not found in the "
            << fvMeshMovers::multiValveEngine::typeName << " mover for region "
            << fluidRegionName_ << nl << nl
            << "Valid valve names are:" << nl << valveNames
            << exit(FatalError);
    }

    return
        septernion
        (
            fluidEngineMover.valves[valveIndex_].isOpen()
          ? fluidEngineMover.valves[valveIndex_].lift()
           *fluidEngineMover.valves[valveIndex_].axis
          : vector::zero
        );
}


bool Foam::solidBodyMotionFunctions::multiValveEngineValveMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    fluidRegionName_ = SBMFCoeffs.lookup<word>("fluidRegion");

    valveName_ = SBMFCoeffs.lookup<word>("valve");

    valveIndex_ = -1;

    return true;
}


// ************************************************************************* //
