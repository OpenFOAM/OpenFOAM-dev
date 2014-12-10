/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "patchToPoly2DMesh.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToPoly2DMesh::flipFaceOrder()
{
    const edgeList& edges = patch_.edges();
    const faceList& localFaces = patch_.localFaces();
    const labelList& meshPoints = patch_.meshPoints();

    Info<< "Flipping face order if necessary." << endl;
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        faces_[edgeI].setSize(2);

        label edgeOwner = owner_[edgeI];

        const face& f = localFaces[edgeOwner];

        label fp = findIndex(f, e[0]);

        if (f.nextLabel(fp) != e[1])
        {
            Info<< "Flipping face " << faces_[edgeI] << endl;
            faces_[edgeI][0] = meshPoints[e[1]];
            faces_[edgeI][1] = meshPoints[e[0]];
        }
        else
        {
            faces_[edgeI][0] = meshPoints[e[0]];
            faces_[edgeI][1] = meshPoints[e[1]];
        }
    }
}


void Foam::patchToPoly2DMesh::createNeighbours()
{
    const edgeList& edges = patch_.edges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    Info<< "Calculating neighbours." << endl;
    forAll(edges, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];
        if (eFaces.size() == 2)
        {
            if (owner_[edgeI] == eFaces[0])
            {
                neighbour_[edgeI] = eFaces[1];
            }
            else
            {
                neighbour_[edgeI] = eFaces[0];
            }
        }
        else if (eFaces.size() == 1)
        {
            continue;
        }
        else
        {
            FatalErrorIn("polyMesh neighbour construction")
                << abort(FatalError);
        }
    }
}


Foam::labelList Foam::patchToPoly2DMesh::internalFaceOrder()
{
    const labelListList& faceEdges = patch_.faceEdges();

    labelList oldToNew(owner_.size(), -1);

    label newFaceI = 0;

    forAll(faceEdges, faceI)
    {
        const labelList& fEdges = faceEdges[faceI];
        // Neighbouring faces
        SortableList<label> nbr(fEdges.size(), -1);

        forAll(fEdges, feI)
        {
            if (fEdges[feI] < neighbour_.size())
            {
                // Internal edge. Get the face on other side.

                label nbrFaceI = neighbour_[fEdges[feI]];

                if (nbrFaceI == faceI)
                {
                    nbrFaceI = owner_[fEdges[feI]];
                }

                if (faceI < nbrFaceI)
                {
                    // faceI is master
                    nbr[feI] = nbrFaceI;
                }
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNew[fEdges[nbr.indices()[i]]] = newFaceI++;
            }
        }
    }

    return oldToNew;
}


void Foam::patchToPoly2DMesh::addPatchFacesToFaces()
{
    const labelList& meshPoints = patch_.meshPoints();

    label offset = patch_.nInternalEdges();
    face f(2);

    forAll(patchNames_, patchI)
    {
        forAllConstIter(EdgeMap<label>, mapEdgesRegion_, eIter)
        {
            if (eIter() == patchI)
            {
                f[0] = meshPoints[eIter.key().start()];
                f[1] = meshPoints[eIter.key().end()];
                faces_[offset++] = f;
            }
        }
    }

    f.clear();
}


void Foam::patchToPoly2DMesh::addPatchFacesToOwner()
{
    const label nInternalEdges = patch_.nInternalEdges();
    const faceList& faces = patch_.faces();
    const label nExternalEdges = patch_.edges().size() - nInternalEdges;
    const labelList& meshPoints = patch_.meshPoints();

    // Reorder patch faces on owner list.
    labelList newOwner = owner_;

    label nMatched = 0;

    for
    (
        label bFaceI = nInternalEdges;
        bFaceI < faces_.size();
        ++bFaceI
    )
    {
        const face& e = faces_[bFaceI];

        bool matched = false;

        for
        (
            label bEdgeI = nInternalEdges;
            bEdgeI < faces_.size();
            ++bEdgeI
        )
        {
            if
            (
                e[0] == meshPoints[patch_.edges()[bEdgeI][0]]
             && e[1] == meshPoints[patch_.edges()[bEdgeI][1]]
            )
            {
                const face& f = faces[owner_[bEdgeI]];

                label fp = findIndex(f, e[0]);

                newOwner[bFaceI] = owner_[bEdgeI];

                if (f.nextLabel(fp) != e[1])
                {
                    Info<< "Flipping" << endl;

                    faces_[bFaceI][0] = e[1];
                    faces_[bFaceI][1] = e[0];
                }

                nMatched++;

                matched = true;
            }
            else if
            (
                e[0] == meshPoints[patch_.edges()[bEdgeI][1]]
             && e[1] == meshPoints[patch_.edges()[bEdgeI][0]]
            )
            {
                Info<< "Warning: Wrong orientation." << endl;
                nMatched++;
                matched = true;
            }
        }
        if (!matched)
        {
            Info<< "No match for edge." << endl;
        }
    }

    if (nMatched != nExternalEdges)
    {
        Info<< "Number of matched edges, " << nMatched
            << ", does not match number of external edges, "
            << nExternalEdges << endl;
    }

    owner_ = newOwner.xfer();
}


void Foam::patchToPoly2DMesh::createPolyMeshComponents()
{
    flipFaceOrder();

    createNeighbours();

    // New function for returning a map of old faces to new faces.
    labelList oldToNew = internalFaceOrder();

    inplaceReorder(oldToNew, faces_);
    inplaceReorder(oldToNew, owner_);
    inplaceReorder(oldToNew, neighbour_);

    // Add patches.
    addPatchFacesToFaces();

    addPatchFacesToOwner();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPoly2DMesh::patchToPoly2DMesh
(
    const MeshedSurface<face>& patch,
    const wordList& patchNames,
    const labelList& patchSizes,
    const EdgeMap<label>& mapEdgesRegion
)
:
    patch_(patch),
    patchNames_(patchNames),
    patchSizes_(patchSizes),
    patchStarts_(patchNames.size(), 0),
    mapEdgesRegion_(mapEdgesRegion),
    points_(patch.points()),
    faces_(patch.nEdges()),
    owner_(PatchTools::edgeOwner(patch)),
    neighbour_(patch.nInternalEdges())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPoly2DMesh::~patchToPoly2DMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::patchToPoly2DMesh::createMesh()
{
    for (label edgeI = 0; edgeI < patch_.nInternalEdges(); edgeI++)
    {
        if (patch_.edgeFaces()[edgeI].size() != 2)
        {
            FatalErrorIn("patchToPoly2DMesh::patchToPoly2DMesh(..)")
                << "internal edge:" << edgeI
                << " patch.edgeFaces()[edgeI]:" << patch_.edgeFaces()[edgeI]
                << abort(FatalError);
        }
    }

    for
    (
        label edgeI = patch_.nInternalEdges();
        edgeI < patch_.nEdges();
        edgeI++
    )
    {
        if (patch_.edgeFaces()[edgeI].size() != 1)
        {
            FatalErrorIn("patchToPoly2DMesh::patchToPoly2DMesh(..)")
                << "boundary edge:" << edgeI
                << " patch.edgeFaces()[edgeI]:" << patch_.edgeFaces()[edgeI]
                << abort(FatalError);
        }
    }

    createPolyMeshComponents();

    label startFace = patch_.nInternalEdges();
    forAll(patchNames_, patchI)
    {
        patchStarts_[patchI] = startFace;
        startFace += patchSizes_[patchI];
    }
}


// ************************************************************************* //
