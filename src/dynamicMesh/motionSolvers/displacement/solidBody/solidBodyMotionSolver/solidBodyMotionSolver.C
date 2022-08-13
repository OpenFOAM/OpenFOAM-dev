/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    pointIDs_(),
    moveAllCells_(false),
    transform_(SBMFPtr_().transformation())
{
    const labelList cellIDs(polyCellSet(mesh, dict).cells());

    moveAllCells_ =
        returnReduce(cellIDs.size() == mesh.nCells(), andOp<bool>());

    if (moveAllCells_)
    {
        Info<< "Applying solid body motion to entire mesh" << endl;
    }
    else
    {
        // collect point IDs of points in cell zone

        boolList movePts(mesh.nPoints(), false);

        forAll(cellIDs, i)
        {
            const cell& c = mesh.cells()[cellIDs[i]];
            forAll(c, j)
            {
                const face& f = mesh.faces()[c[j]];
                forAll(f, k)
                {
                    movePts[f[k]] = true;
                }
            }
        }

        syncTools::syncPointList(mesh, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs(mesh.nPoints());
        forAll(movePts, i)
        {
            if (movePts[i])
            {
                ptIDs.append(i);
            }
        }

        pointIDs_.transfer(ptIDs);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionSolver::~solidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::solidBodyMotionSolver::curPoints() const
{
    transform_ = SBMFPtr_().transformation();

    if (moveAllCells_)
    {
        return transformPoints(transform_, points0_);
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs_) = transformPoints
        (
            transform_,
            pointField(points0_, pointIDs_)
        );

        return ttransformedPts;
    }
}


void Foam::solidBodyMotionSolver::topoChange(const polyTopoChangeMap& map)
{
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
        label oldPointi = map.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            const label masterPointi = map.reversePointMap()[oldPointi];

            if (masterPointi == pointi)
            {
                newPoints0[pointi] = points0_[oldPointi];
            }
            else
            {
                newPoints0[pointi] =
                    transform_.invTransformPoint(points[pointi]);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot determine co-ordinates of introduced vertices."
                << " New vertex " << pointi << " at co-ordinate "
                << points[pointi] << exit(FatalError);
        }
    }

    twoDCorrectPoints(newPoints0);

    points0_.transfer(newPoints0);
}


// ************************************************************************* //
