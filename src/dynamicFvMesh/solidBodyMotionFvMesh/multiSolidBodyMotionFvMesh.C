/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "multiSolidBodyMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "cellZoneMesh.H"
#include "boolList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        multiSolidBodyMotionFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionFvMesh::multiSolidBodyMotionFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    )
{
    if (undisplacedPoints_.size() != nPoints())
    {
        FatalIOErrorIn
        (
            "multiSolidBodyMotionFvMesh::multiSolidBodyMotionFvMesh"
            "(const IOobject&)",
            dynamicMeshCoeffs_
        )   << "Read " << undisplacedPoints_.size()
            << " undisplaced points from " << undisplacedPoints_.objectPath()
            << " but the current mesh has " << nPoints()
            << exit(FatalIOError);
    }


    zoneIDs_.setSize(dynamicMeshCoeffs_.size());
    SBMFs_.setSize(dynamicMeshCoeffs_.size());
    pointIDs_.setSize(dynamicMeshCoeffs_.size());
    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs_, iter)
    {
        if (iter().isDict())
        {
            zoneIDs_[zoneI] = cellZones().findZoneID(iter().keyword());

            if (zoneIDs_[zoneI] == -1)
            {
                FatalIOErrorIn
                (
                    "multiSolidBodyMotionFvMesh::"
                    "multiSolidBodyMotionFvMesh(const IOobject&)",
                    dynamicMeshCoeffs_
                )   << "Cannot find cellZone named " << iter().keyword()
                    << ". Valid zones are " << cellZones().names()
                    << exit(FatalIOError);
            }

            const dictionary& subDict = iter().dict();

            SBMFs_.set
            (
                zoneI,
                solidBodyMotionFunction::New(subDict, io.time())
            );

            // Collect points of cell zone.
            const cellZone& cz = cellZones()[zoneIDs_[zoneI]];

            boolList movePts(nPoints(), false);

            forAll(cz, i)
            {
                label cellI = cz[i];
                const cell& c = cells()[cellI];
                forAll(c, j)
                {
                    const face& f = faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointI = f[k];
                        movePts[pointI] = true;
                    }
                }
            }

            syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zoneI].transfer(ptIDs);

            Info<< "Applying solid body motion " << SBMFs_[zoneI].type()
                << " to " << pointIDs_[zoneI].size() << " points of cellZone "
                << iter().keyword() << endl;

            zoneI++;
        }
    }
    zoneIDs_.setSize(zoneI);
    SBMFs_.setSize(zoneI);
    pointIDs_.setSize(zoneI);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionFvMesh::~multiSolidBodyMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiSolidBodyMotionFvMesh::update()
{
    static bool hasWarned = false;

    pointField transformedPts(undisplacedPoints_);

    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];

        UIndirectList<point>(transformedPts, zonePoints) =
            transform
            (
                SBMFs_[i].transformation(),
                pointField(transformedPts, zonePoints)
            );
    }

    fvMesh::movePoints(transformedPts);

    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningIn("multiSolidBodyMotionFvMesh::update()")
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //
