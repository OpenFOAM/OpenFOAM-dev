/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
#include "cellZoneMesh.H"
#include "cellSet.H"
#include "boolList.H"
#include "syncTools.H"

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
    moveAllCells_(false)
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
    if (moveAllCells_)
    {
        return transformPoints(SBMFPtr_().transformation(), points0_);
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs_) = transformPoints
        (
            SBMFPtr_().transformation(),
            pointField(points0_, pointIDs_)
        );

        return ttransformedPts;
    }
}


// ************************************************************************* //
