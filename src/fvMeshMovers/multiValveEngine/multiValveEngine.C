/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "multiValveEngine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(multiValveEngine, 0);
    addToRunTimeSelectionTable(fvMeshMover, multiValveEngine, fvMesh);
}
}

Foam::word Foam::fvMeshMovers::multiValveEngine::cylinderHeadName
(
    "cylinderHead"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::labelHashSet
Foam::fvMeshMovers::multiValveEngine::findLinerPatchSet() const
{
    return mesh().boundaryMesh().patchSet
    (
        wordReList(dict().lookup("linerPatches"))
    );
}


Foam::labelHashSet Foam::fvMeshMovers::multiValveEngine::findSlidingPatchSet()
{
    return mesh().boundaryMesh().patchSet
    (
        wordReList(dict().lookup("slidingPatches"))
    );
}


Foam::labelHashSet Foam::fvMeshMovers::multiValveEngine::findStaticPatchSet()
{
    labelHashSet movingPatchSet(piston_.patchSet_);

    forAll(valves_, valvei)
    {
        movingPatchSet += valves_[valvei].patchSet_;
    }

    labelHashSet staticPatchSet;

    forAll(mesh().boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        // Exclude non-static patches
        if
        (
           !polyPatch::constraintType(pp.type())
        && !slidingPatchSet_.found(pp.index())
        && !movingPatchSet.found(pp.index())
        )
        {
            staticPatchSet.insert(pp.index());
        }
    }

    return staticPatchSet;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiValveEngine::multiValveEngine(fvMesh& mesh)
:
    fvMeshMover(mesh),
    linerPatchSet_(findLinerPatchSet()),
    slidingPatchSet_(findSlidingPatchSet()),
    piston_("piston", *this, dict().subDict("piston")),
    valves_(*this, dict().subOrEmptyDict("valves")),
    staticPatchSet_(findStaticPatchSet()),
    frozenPointZones_
    (
        dict().lookupOrDefault("frozenZones", wordReList::null())
    ),
    linerPatchSet(linerPatchSet_),
    slidingPatchSet(slidingPatchSet_),
    piston(piston_),
    valves(valves_),
    staticPatchSet(staticPatchSet_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiValveEngine::~multiValveEngine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fvMeshMovers::multiValveEngine::userTime() const
{
    return mesh().time().userTimeValue();
}


Foam::scalar Foam::fvMeshMovers::multiValveEngine::userDeltaT() const
{
    return mesh().time().userDeltaTValue();
}


bool Foam::fvMeshMovers::multiValveEngine::update()
{
    // Accumulate point motion from each moving object into newPoints field.
    pointField newPoints(mesh().points());

    piston_.updatePoints(newPoints);

    forAll(valves_, valvei)
    {
        valves_[valvei].updatePoints(newPoints);
    }

    // Update the mesh according to the newPoints field.
    mesh().movePoints(newPoints);

    return true;
}


void Foam::fvMeshMovers::multiValveEngine::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fvMeshMovers::multiValveEngine::mapMesh(const polyMeshMap& map)
{
    slidingPatchSet_ = findSlidingPatchSet();
    linerPatchSet_ = findLinerPatchSet();
    staticPatchSet_ = findStaticPatchSet();

    piston_.mapMesh(map);

    forAll(valves_, valvei)
    {
        valves_[valvei].mapMesh(map);
    }
}


void Foam::fvMeshMovers::multiValveEngine::distribute
(
    const polyDistributionMap&
)
{
    NotImplemented;
}


// ************************************************************************* //
