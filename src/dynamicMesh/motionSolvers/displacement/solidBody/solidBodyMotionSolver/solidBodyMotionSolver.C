/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "meshCellZones.H"
#include "cellSet.H"
#include "boolList.H"
#include "syncTools.H"
#include "mapPolyMesh.H"

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
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(mesh, dict, typeName),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    pointIDs_(),
    moveAllCells_(false),
    transform_(SBMFPtr_().transformation())
{
    word cellZoneName =
        coeffDict().lookupOrDefault<word>("cellZone", "none");

    word cellSetName =
        coeffDict().lookupOrDefault<word>("cellSet", "none");

    if ((cellZoneName != "none") && (cellSetName != "none"))
    {
        FatalIOErrorInFunction(coeffDict())
            << "Either cellZone OR cellSet can be supplied, but not both. "
            << "If neither is supplied, all cells will be included"
            << exit(FatalIOError);
    }

    labelList cellIDs;
    if (cellZoneName != "none")
    {
        Info<< "Applying solid body motion to cellZone " << cellZoneName
            << endl;

        label zoneID = mesh.cellZones().findZoneID(cellZoneName);

        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Unable to find cellZone " << cellZoneName
                << ".  Valid cellZones are:"
                << mesh.cellZones().names()
                << exit(FatalError);
        }

        cellIDs = mesh.cellZones()[zoneID];
    }

    if (cellSetName != "none")
    {
        Info<< "Applying solid body motion to cellSet " << cellSetName
            << endl;

        cellSet set(mesh, cellSetName);

        cellIDs = set.toc();
    }

    label nCells = returnReduce(cellIDs.size(), sumOp<label>());
    moveAllCells_ = nCells == 0;

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
            label celli = cellIDs[i];
            const cell& c = mesh.cells()[celli];
            forAll(c, j)
            {
                const face& f = mesh.faces()[c[j]];
                forAll(f, k)
                {
                    label pointi = f[k];
                    movePts[pointi] = true;
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


void Foam::solidBodyMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // pointMesh already updates pointFields

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        mpm.hasMotionPoints()
      ? mpm.preMotionPoints()
      : mesh().points()
    );

    pointField newPoints0(mpm.pointMap().size());

    forAll(newPoints0, pointi)
    {
        label oldPointi = mpm.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            const label masterPointi = mpm.reversePointMap()[oldPointi];

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

    // points0 changed - set to write and check-in to database
    points0_.rename("points0");
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().timeName();
    points0_.checkIn();
}


// ************************************************************************* //
