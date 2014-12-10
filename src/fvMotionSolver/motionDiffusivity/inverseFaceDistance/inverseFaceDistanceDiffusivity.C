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

#include "inverseFaceDistanceDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"
#include "wallPoint.H"
#include "MeshWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseFaceDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseFaceDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::inverseFaceDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    patchNames_(mdData)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::~inverseFaceDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inverseFaceDistanceDiffusivity::correct()
{
    const polyBoundaryMesh& bdry = mesh().boundaryMesh();

    labelHashSet patchSet(bdry.size());

    label nPatchFaces = 0;

    forAll(patchNames_, i)
    {
        const label pID = bdry.findPatchID(patchNames_[i]);

        if (pID > -1)
        {
            patchSet.insert(pID);
            nPatchFaces += bdry[pID].size();
        }
    }

    List<wallPoint> faceDist(nPatchFaces);
    labelList changedFaces(nPatchFaces);

    nPatchFaces = 0;

    forAllConstIter(labelHashSet, patchSet, iter)
    {
        const polyPatch& patch = bdry[iter.key()];

        const vectorField::subField fc(patch.faceCentres());

        forAll(fc, patchFaceI)
        {
            changedFaces[nPatchFaces] = patch.start() + patchFaceI;

            faceDist[nPatchFaces] = wallPoint(fc[patchFaceI], 0);

            nPatchFaces++;
        }
    }
    faceDist.setSize(nPatchFaces);
    changedFaces.setSize(nPatchFaces);

    MeshWave<wallPoint> waveInfo
    (
        mesh(),
        changedFaces,
        faceDist,
        mesh().globalData().nTotalCells()+1   // max iterations
    );

    const List<wallPoint>& faceInfo = waveInfo.allFaceInfo();
    const List<wallPoint>& cellInfo = waveInfo.allCellInfo();

    for (label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        scalar dist = faceInfo[faceI].distSqr();

        faceDiffusivity_[faceI] = 1.0/sqrt(dist);
    }

    forAll(faceDiffusivity_.boundaryField(), patchI)
    {
        fvsPatchScalarField& bfld = faceDiffusivity_.boundaryField()[patchI];

        const labelUList& faceCells = bfld.patch().faceCells();

        if (patchSet.found(patchI))
        {
            forAll(bfld, i)
            {
                scalar dist = cellInfo[faceCells[i]].distSqr();
                bfld[i] = 1.0/sqrt(dist);
            }
        }
        else
        {
            const label start = bfld.patch().start();

            forAll(bfld, i)
            {
                scalar dist = faceInfo[start+i].distSqr();
                bfld[i] = 1.0/sqrt(dist);
            }
        }
    }
}


// ************************************************************************* //
