/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "multiSolidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "cellZoneList.H"
#include "boolList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        multiSolidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::multiSolidBodyMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName)
{
    zoneIndices_.setSize(dict.size());
    SBMFs_.setSize(dict.size());
    pointIndices_.setSize(dict.size());
    label zonei = 0;

    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            zoneIndices_[zonei] = mesh.cellZones().findIndex(iter().keyword());

            if (zoneIndices_[zonei] == -1)
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Cannot find cellZone named " << iter().keyword()
                    << ". Valid zones are " << mesh.cellZones().toc()
                    << exit(FatalIOError);
            }

            const dictionary& subDict = iter().dict();

            SBMFs_.set
            (
                zonei,
                solidBodyMotionFunction::New(subDict, mesh.time())
            );

            // Collect points of cell zone.
            const cellZone& cz = mesh.cellZones()[zoneIndices_[zonei]];

            boolList movePts(mesh.nPoints(), false);

            forAll(cz, i)
            {
                label celli = cz[i];
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

            pointIndices_[zonei].transfer(ptIDs);

            Info<< "Applying solid body motion " << SBMFs_[zonei].type()
                << " to " << pointIndices_[zonei].size()
                << " points of cellZone " << iter().keyword()
                << endl;

            zonei++;
        }
    }
    zoneIndices_.setSize(zonei);
    SBMFs_.setSize(zonei);
    pointIndices_.setSize(zonei);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::~multiSolidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::multiSolidBodyMotionSolver::curPoints() const
{
    tmp<pointField> ttransformedPts(new pointField(mesh().points()));
    pointField& transformedPts = ttransformedPts.ref();

    forAll(zoneIndices_, i)
    {
        const labelList& zonePoints = pointIndices_[i];

        UIndirectList<point>(transformedPts, zonePoints) = transformPoints
        (
            SBMFs_[i].transformation(),
            pointField(points0_, zonePoints)
        );
    }

    return ttransformedPts;
}


// ************************************************************************* //
