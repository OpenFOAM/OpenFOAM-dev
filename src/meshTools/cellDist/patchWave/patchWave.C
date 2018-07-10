/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "patchWave.H"
#include "polyMesh.H"
#include "wallPoint.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchWave::setChangedFaces
(
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<wallPoint>& faceDist
) const
{
    const polyMesh& mesh = cellDistFuncs::mesh();

    label nChangedFaces = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            forAll(patch.faceCentres(), patchFacei)
            {
                label meshFacei = patch.start() + patchFacei;

                changedFaces[nChangedFaces] = meshFacei;

                faceDist[nChangedFaces] =
                    wallPoint
                    (
                        patch.faceCentres()[patchFacei],
                        0.0
                    );

                nChangedFaces++;
            }
        }
    }
}


Foam::label Foam::patchWave::getValues(const MeshWave<wallPoint>& waveInfo)
{
    const List<wallPoint>& cellInfo = waveInfo.allCellInfo();
    const List<wallPoint>& faceInfo = waveInfo.allFaceInfo();

    label nIllegal = 0;

    // Copy cell values
    distance_.setSize(cellInfo.size());

    forAll(cellInfo, celli)
    {
        scalar dist = cellInfo[celli].distSqr();

        if (cellInfo[celli].valid(waveInfo.data()))
        {
            distance_[celli] = Foam::sqrt(dist);
        }
        else
        {
            distance_[celli] = dist;

            nIllegal++;
        }
    }

    // Copy boundary values
    forAll(patchDistance_, patchi)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchi];

        // Allocate storage for patchDistance
        scalarField* patchDistPtr = new scalarField(patch.size());

        patchDistance_.set(patchi, patchDistPtr);

        scalarField& patchField = *patchDistPtr;

        forAll(patchField, patchFacei)
        {
            label meshFacei = patch.start() + patchFacei;

            scalar dist = faceInfo[meshFacei].distSqr();

            if (faceInfo[meshFacei].valid(waveInfo.data()))
            {
                // Adding small to avoid problems with /0 in the turbulence
                // models
                patchField[patchFacei] = Foam::sqrt(dist) + small;
            }
            else
            {
                patchField[patchFacei] = dist;

                nIllegal++;
            }
        }
    }
    return nIllegal;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchWave::patchWave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    const bool correctWalls
)
:
    cellDistFuncs(mesh),
    patchIDs_(patchIDs),
    correctWalls_(correctWalls),
    nUnset_(0),
    distance_(mesh.nCells()),
    patchDistance_(mesh.boundaryMesh().size())
{
    patchWave::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchWave::~patchWave()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchWave::correct()
{
    // Set initial changed faces: set wallPoint for wall faces to wall centre

    label nPatch = sumPatchSize(patchIDs_);

    List<wallPoint> faceDist(nPatch);
    labelList changedFaces(nPatch);

    // Set to faceDist information to facecentre on walls.
    setChangedFaces(patchIDs_, changedFaces, faceDist);

    // Do calculate wall distance by 'growing' from faces.
    MeshWave<wallPoint> waveInfo
    (
        mesh(),
        changedFaces,
        faceDist,
        mesh().globalData().nTotalCells()+1 // max iterations
    );

    // Copy distance into return field
    nUnset_ = getValues(waveInfo);

    // Correct wall cells for true distance
    if (correctWalls_)
    {
        Map<label> nearestFace(2*nPatch);

        correctBoundaryFaceCells
        (
            patchIDs_,
            distance_,
            nearestFace
        );

        correctBoundaryPointCells
        (
            patchIDs_,
            distance_,
            nearestFace
        );
    }
}


// ************************************************************************* //
