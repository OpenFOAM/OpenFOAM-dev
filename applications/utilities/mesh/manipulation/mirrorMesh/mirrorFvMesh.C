/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "mirrorFvMesh.H"
#include "Time.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mirrorFvMesh::mirrorFvMesh(const IOobject& io)
:
    fvMesh(io),
    mirrorMeshDict_
    (
        IOobject
        (
            "mirrorMeshDict",
            time().system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mirrorMeshPtr_(NULL)
{
    plane mirrorPlane(mirrorMeshDict_);

    scalar planeTolerance
    (
        readScalar(mirrorMeshDict_.lookup("planeTolerance"))
    );

    const pointField& oldPoints = points();
    const faceList& oldFaces = faces();
    const cellList& oldCells = cells();
    const label nOldInternalFaces = nInternalFaces();
    const polyPatchList& oldPatches = boundaryMesh();

    // Mirror the points
    Info<< "Mirroring points. Old points: " << oldPoints.size();

    pointField newPoints(2*oldPoints.size());
    label nNewPoints = 0;

    labelList mirrorPointLookup(oldPoints.size(), -1);

    // Grab the old points
    forAll(oldPoints, pointI)
    {
        newPoints[nNewPoints] = oldPoints[pointI];
        nNewPoints++;
    }

    forAll(oldPoints, pointI)
    {
        scalar alpha =
            mirrorPlane.normalIntersect
            (
                oldPoints[pointI],
                mirrorPlane.normal()
            );

        // Check plane on tolerance
        if (mag(alpha) > planeTolerance)
        {
            // The point gets mirrored
            newPoints[nNewPoints] =
                oldPoints[pointI] + 2.0*alpha*mirrorPlane.normal();

            // remember the point correspondence
            mirrorPointLookup[pointI] = nNewPoints;
            nNewPoints++;
        }
        else
        {
            // The point is on the plane and does not get mirrored
            // Adjust plane location
            newPoints[nNewPoints] =
                oldPoints[pointI] + alpha*mirrorPlane.normal();

            mirrorPointLookup[pointI] = pointI;
        }
    }

    // Reset the size of the point list
    Info<< " New points: " << nNewPoints << endl;
    newPoints.setSize(nNewPoints);

    Info<< "Mirroring faces. Old faces: " << oldFaces.size();

    // Algorithm:
    // During mirroring, the faces that were previously boundary faces
    // in the mirror plane may become ineternal faces. In order to
    // deal with the ordering of the faces, the algorithm is split
    // into two parts.  For original faces, the internal faces are
    // distributed to their owner cells.  Once all internal faces are
    // distributed, the boundary faces are visited and if they are in
    // the mirror plane they are added to the master cells (the future
    // boundary faces are not touched).  After the first phase, the
    // internal faces are collected in the cell order and numbering
    // information is added.  Then, the internal faces are mirrored
    // and the face numbering data is stored for the mirrored section.
    // Once all the internal faces are mirrored, the boundary faces
    // are added by mirroring the faces patch by patch.

    // Distribute internal faces
    labelListList newCellFaces(oldCells.size());

    const labelUList& oldOwnerStart = lduAddr().ownerStartAddr();

    forAll(newCellFaces, cellI)
    {
        labelList& curFaces = newCellFaces[cellI];

        const label s = oldOwnerStart[cellI];
        const label e = oldOwnerStart[cellI + 1];

        curFaces.setSize(e - s);

        forAll(curFaces, i)
        {
            curFaces[i] = s + i;
        }
    }

    // Distribute boundary faces.  Remember the faces that have been inserted
    // as internal
    boolListList insertedBouFace(oldPatches.size());

    forAll(oldPatches, patchI)
    {
        const polyPatch& curPatch = oldPatches[patchI];

        if (curPatch.coupled())
        {
            WarningIn("mirrorFvMesh::mirrorFvMesh(const IOobject&)")
                << "Found coupled patch " << curPatch.name() << endl
                << "    Mirroring faces on coupled patches destroys"
                << " the ordering. This might be fixed by running a dummy"
                << " createPatch afterwards." << endl;
        }

        boolList& curInsBouFace = insertedBouFace[patchI];

        curInsBouFace.setSize(curPatch.size());
        curInsBouFace = false;

        // Get faceCells for face insertion
        const labelUList& curFaceCells = curPatch.faceCells();
        const label curStart = curPatch.start();

        forAll(curPatch, faceI)
        {
            // Find out if the mirrored face is identical to the
            // original.  If so, the face needs to become internal and
            // added to its owner cell
            const face& origFace = curPatch[faceI];

            face mirrorFace(origFace.size());
            forAll(mirrorFace, pointI)
            {
                mirrorFace[pointI] = mirrorPointLookup[origFace[pointI]];
            }

            if (origFace == mirrorFace)
            {
                // The mirror is identical to current face.  This will
                // become an internal face
                const label oldSize = newCellFaces[curFaceCells[faceI]].size();

                newCellFaces[curFaceCells[faceI]].setSize(oldSize + 1);
                newCellFaces[curFaceCells[faceI]][oldSize] = curStart + faceI;

                curInsBouFace[faceI] = true;
            }
        }
    }

    // Construct the new list of faces.  Boundary faces are added
    // last, cush that each patch is mirrored separately.  The
    // addressing is stored in two separate arrays: first for the
    // original cells (face order has changed) and then for the
    // mirrored cells.
    labelList masterFaceLookup(oldFaces.size(), -1);
    labelList mirrorFaceLookup(oldFaces.size(), -1);

    faceList newFaces(2*oldFaces.size());
    label nNewFaces = 0;

    // Insert original (internal) faces
    forAll(newCellFaces, cellI)
    {
        const labelList& curCellFaces = newCellFaces[cellI];

        forAll(curCellFaces, cfI)
        {
            newFaces[nNewFaces] = oldFaces[curCellFaces[cfI]];
            masterFaceLookup[curCellFaces[cfI]] = nNewFaces;

            nNewFaces++;
        }
    }

    // Mirror internal faces
    for (label faceI = 0; faceI < nOldInternalFaces; faceI++)
    {
        const face& oldFace = oldFaces[faceI];
        face& nf = newFaces[nNewFaces];
        nf.setSize(oldFace.size());

        nf[0] = mirrorPointLookup[oldFace[0]];

        for (label i = 1; i < oldFace.size(); i++)
        {
            nf[i] = mirrorPointLookup[oldFace[oldFace.size() - i]];
        }

        mirrorFaceLookup[faceI] = nNewFaces;
        nNewFaces++;
    }

    // Mirror boundary faces patch by patch

    wordList newPatchTypes(boundary().size());
    wordList newPatchNames(boundary().size());
    labelList newPatchSizes(boundary().size(), -1);
    labelList newPatchStarts(boundary().size(), -1);
    label nNewPatches = 0;

    forAll(boundaryMesh(), patchI)
    {
        const label curPatchSize = boundaryMesh()[patchI].size();
        const label curPatchStart = boundaryMesh()[patchI].start();
        const boolList& curInserted = insertedBouFace[patchI];

        newPatchStarts[nNewPatches] = nNewFaces;

        // Master side
        for (label faceI = 0; faceI < curPatchSize; faceI++)
        {
            // Check if the face has already been added.  If not, add it and
            // insert the numbering details.
            if (!curInserted[faceI])
            {
                newFaces[nNewFaces] = oldFaces[curPatchStart + faceI];

                masterFaceLookup[curPatchStart + faceI] = nNewFaces;
                nNewFaces++;
            }
        }

        // Mirror side
        for (label faceI = 0; faceI < curPatchSize; faceI++)
        {
            // Check if the face has already been added.  If not, add it and
            // insert the numbering details.
            if (!curInserted[faceI])
            {
                const face& oldFace = oldFaces[curPatchStart + faceI];
                face& nf = newFaces[nNewFaces];
                nf.setSize(oldFace.size());

                nf[0] = mirrorPointLookup[oldFace[0]];

                for (label i = 1; i < oldFace.size(); i++)
                {
                    nf[i] = mirrorPointLookup[oldFace[oldFace.size() - i]];
                }

                mirrorFaceLookup[curPatchStart + faceI] = nNewFaces;
                nNewFaces++;
            }
            else
            {
                // Grab the index of the master face for the mirror side
                mirrorFaceLookup[curPatchStart + faceI] =
                    masterFaceLookup[curPatchStart + faceI];
            }
        }

        // If patch exists, grab the name and type of the original patch
        if (nNewFaces > newPatchStarts[nNewPatches])
        {
            newPatchTypes[nNewPatches] = boundaryMesh()[patchI].type();
            newPatchNames[nNewPatches] = boundaryMesh()[patchI].name();
            newPatchSizes[nNewPatches] =
                nNewFaces - newPatchStarts[nNewPatches];

            nNewPatches++;
        }
    }

    // Tidy up the lists
    newFaces.setSize(nNewFaces);
    Info<< " New faces: " << nNewFaces << endl;

    newPatchTypes.setSize(nNewPatches);
    newPatchNames.setSize(nNewPatches);
    newPatchSizes.setSize(nNewPatches);
    newPatchStarts.setSize(nNewPatches);

    Info<< "Mirroring patches. Old patches: " << boundary().size()
        << " New patches: " << nNewPatches << endl;

    Info<< "Mirroring cells.  Old cells: " << oldCells.size()
        << " New cells: " << 2*oldCells.size() << endl;

    cellList newCells(2*oldCells.size());
    label nNewCells = 0;

    // Grab the original cells.  Take care of face renumbering.
    forAll(oldCells, cellI)
    {
        const cell& oc = oldCells[cellI];

        cell& nc = newCells[nNewCells];
        nc.setSize(oc.size());

        forAll(oc, i)
        {
            nc[i] = masterFaceLookup[oc[i]];
        }

        nNewCells++;
    }

    // Mirror the cells
    forAll(oldCells, cellI)
    {
        const cell& oc = oldCells[cellI];

        cell& nc = newCells[nNewCells];
        nc.setSize(oc.size());

        forAll(oc, i)
        {
            nc[i] = mirrorFaceLookup[oc[i]];
        }

        nNewCells++;
    }

    // Mirror the cell shapes
    Info<< "Mirroring cell shapes." << endl;

    Info<< nl << "Creating new mesh" << endl;
    mirrorMeshPtr_ = new fvMesh
    (
        io,
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newCells)
    );

    fvMesh& pMesh = *mirrorMeshPtr_;

    // Add the boundary patches
    List<polyPatch*> p(newPatchTypes.size());

    forAll(p, patchI)
    {
        p[patchI] = polyPatch::New
        (
            newPatchTypes[patchI],
            newPatchNames[patchI],
            newPatchSizes[patchI],
            newPatchStarts[patchI],
            patchI,
            pMesh.boundaryMesh()
        ).ptr();
    }

    pMesh.addPatches(p);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mirrorFvMesh::~mirrorFvMesh()
{}


// ************************************************************************* //
