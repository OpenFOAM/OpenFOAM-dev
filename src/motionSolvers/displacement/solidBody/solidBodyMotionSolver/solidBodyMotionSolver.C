/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "solidBodyMotionSolver.H"
#include "polyCellSet.H"
#include "transformField.H"
#include "syncTools.H"
#include "polyTopoChangeMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        solidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidBodyMotionSolver::updateSetPointIDs()
{
    if (set_.selectionType() == polyCellSet::selectionTypes::all)
    {
        setPointIDs_.clear();
        return;
    }

    boolList pointInSet(mesh().nPoints(), false);

    forAll(set_.cells(), setCelli)
    {
        const cell& c = mesh().cells()[set_.cells()[setCelli]];
        forAll(c, cfi)
        {
            const face& f = mesh().faces()[c[cfi]];
            forAll(f, fpi)
            {
                pointInSet[f[fpi]] = true;
            }
        }
    }

    syncTools::syncPointList(mesh(), pointInSet, orEqOp<bool>(), false);

    setPointIDs_.resize(count(pointInSet, true));

    label setPointi = 0;
    forAll(pointInSet, pointi)
    {
        if (pointInSet[pointi])
        {
            setPointIDs_[setPointi ++] = pointi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionSolver::solidBodyMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    set_(mesh, dict),
    setPointIDs_(),
    transform_(SBMFPtr_().transformation())
{
    if (set_.selectionType() == polyCellSet::selectionTypes::all)
    {
        Info<< "Applying solid body motion to entire mesh" << endl;
    }

    updateSetPointIDs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionSolver::~solidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::solidBodyMotionSolver::curPoints() const
{
    transform_ = SBMFPtr_().transformation();

    if (set_.selectionType() == polyCellSet::selectionTypes::all)
    {
        return transformPoints(transform_, points0_);
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, setPointIDs_) = transformPoints
        (
            transform_,
            pointField(points0_, setPointIDs_)
        );

        return ttransformedPts;
    }
}


void Foam::solidBodyMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
    updateSetPointIDs();

    boolList pointInSet
    (
        mesh().nPoints(),
        set_.selectionType() == polyCellSet::selectionTypes::all
    );
    UIndirectList<bool>(pointInSet, setPointIDs_) = true;

    // pointMesh already updates pointFields

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        map.hasMotionPoints()
      ? map.preMotionPoints()
      : mesh().points()
    );

    pointField newPoints0(map.pointMap().size());

    forAll(newPoints0, pointi)
    {
        const label oldPointi = map.pointMap()[pointi];

        if (oldPointi < 0)
        {
            FatalErrorInFunction
                << "Cannot determine co-ordinates of introduced vertices."
                << " New vertex " << pointi << " at co-ordinate "
                << points[pointi] << exit(FatalError);
        }

        if (map.reversePointMap()[oldPointi] == pointi)
        {
            newPoints0[pointi] = points0_[oldPointi];
        }
        else
        {
            newPoints0[pointi] =
                pointInSet[pointi]
              ? transform_.invTransformPoint(points[pointi])
              : points[pointi];
        }
    }

    twoDCorrectPoints(newPoints0);

    // Move into base class storage and mark as to-be-written
    points0_.transfer(newPoints0);
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().name();
}


void Foam::solidBodyMotionSolver::distribute(const polyDistributionMap& map)
{
    points0MotionSolver::distribute(map);

    set_.distribute(map);
    updateSetPointIDs();
}


void Foam::solidBodyMotionSolver::mapMesh(const polyMeshMap& map)
{
    points0MotionSolver::mapMesh(map);

    set_.mapMesh(map);
    updateSetPointIDs();
}


// ************************************************************************* //
