/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "tetDecomposer.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetDecomposer, 0);

    template<>
    const char* NamedEnum<tetDecomposer::decompositionType, 2>::names[] =
    {
        "faceCentre",
        "faceDiagonal"
    };

    const NamedEnum<tetDecomposer::decompositionType, 2>
        tetDecomposer::decompositionTypeNames;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetDecomposer::modifyFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label faceI,
    const label own,
    const label nei,
    const label patchI,
    const label zoneI,
    const bool zoneFlip
) const
{
    // First usage of face. Modify.
    if (nei == -1 || own < nei)
    {
        meshMod.modifyFace
        (
            f,                          // modified face
            faceI,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patchI,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }
    else
    {
        meshMod.modifyFace
        (
            f.reverseFace(),            // modified face
            faceI,                      // label of face
            nei,                        // owner
            own,                        // neighbour
            true,                       // face flip
            patchI,                     // patch for face
            zoneI,                      // zone for face
            !zoneFlip                   // face flip in zone
        );
    }
}


void Foam::tetDecomposer::addFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label patchI,
    const label zoneI,
    const bool zoneFlip
) const
{
    // Second or more usage of face. Add.
    if (nei == -1 || own < nei)
    {
        meshMod.addFace
        (
            f,                          // modified face
            own,                        // owner
            nei,                        // neighbour
            masterPointID,              // master point
            masterEdgeID,               // master edge
            masterFaceID,               // master face
            false,                      // face flip
            patchI,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }
    else
    {
        meshMod.addFace
        (
            f.reverseFace(),            // modified face
            nei,                        // owner
            own,                        // neighbour
            masterPointID,              // master point
            masterEdgeID,               // master edge
            masterFaceID,               // master face
            true,                       // face flip
            patchI,                     // patch for face
            zoneI,                      // zone for face
            !zoneFlip                   // face flip in zone
        );
    }
}


// Work out triangle index given the starting vertex in the face
Foam::label Foam::tetDecomposer::triIndex(const label faceI, const label fp)
const
{
    const face& f = mesh_.faces()[faceI];
    const label fp0 = mesh_.tetBasePtIs()[faceI];

    // Work out triangle index on this face
    label thisTriI;
    if (fp == fp0)
    {
        thisTriI = 0;
    }
    else if (fp == f.rcIndex(fp0))
    {
        thisTriI = f.size()-3;
    }
    else
    {
        thisTriI = (fp-fp0-1) % (f.size()-2);
    }
    return thisTriI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetDecomposer::tetDecomposer(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetDecomposer::setRefinement
(
    const decompositionType decomposeType,
    polyTopoChange& meshMod
)
{
    cellToPoint_.setSize(mesh_.nCells());
    forAll(mesh_.cellCentres(), cellI)
    {
        // Any point on the cell
        label masterPointI = mesh_.faces()[mesh_.cells()[cellI][0]][0];

        cellToPoint_[cellI] = meshMod.addPoint
        (
            mesh_.cellCentres()[cellI],
            masterPointI,
            -1,
            true
        );
    }


    // Add face centre points
    if (decomposeType == FACECENTRETETS)
    {
        faceToPoint_.setSize(mesh_.nFaces());
        forAll(mesh_.faceCentres(), faceI)
        {
            // Any point on the face
            const label masterPointI = mesh_.faces()[faceI][0];

            faceToPoint_[faceI] = meshMod.addPoint
            (
                mesh_.faceCentres()[faceI],
                masterPointI,
                -1,
                true
            );
        }
    }


    // Per face, per point (faceCentre) or triangle (faceDiag) the added cell
    faceOwnerCells_.setSize(mesh_.nFaces());
    faceNeighbourCells_.setSize(mesh_.nFaces());

    if (decomposeType == FACECENTRETETS)
    {
        forAll(faceOwnerCells_, faceI)
        {
            const face& f = mesh_.faces()[faceI];
            faceOwnerCells_[faceI].setSize(f.size(), -1);
            faceNeighbourCells_[faceI].setSize(f.size(), -1);
        }
    }
    else
    {
        // Force construction of diagonal decomposition
        (void)mesh_.tetBasePtIs();

        forAll(faceOwnerCells_, faceI)
        {
            const face& f = mesh_.faces()[faceI];
            faceOwnerCells_[faceI].setSize(f.size()-2, -1);
            faceNeighbourCells_[faceI].setSize(f.size()-2, -1);
        }
    }


    forAll(mesh_.cells(), cellI)
    {
        const cell& cFaces = mesh_.cells()[cellI];

        EdgeMap<label> edgeToFace(8*cFaces.size());

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const face& f = mesh_.faces()[faceI];

            // Get reference to either owner or neighbour
            labelList& added =
            (
                (mesh_.faceOwner()[faceI] == cellI)
              ? faceOwnerCells_[faceI]
              : faceNeighbourCells_[faceI]
            );

            if (decomposeType == FACECENTRETETS)
            {
                forAll(f, fp)
                {
                    if (cFaceI == 0 && fp == 0)
                    {
                        // Reuse cell itself
                        added[fp] = cellI;
                    }
                    else
                    {
                        added[fp] = meshMod.addCell
                        (
                            -1,     // masterPoint
                            -1,     // masterEdge
                            -1,     // masterFace
                            cellI,  // masterCell
                            mesh_.cellZones().whichZone(cellI)
                        );
                    }
                }
            }
            else
            {
                for (label triI = 0; triI < f.size()-2; triI++)
                {
                    if (cFaceI == 0 && triI == 0)
                    {
                        // Reuse cell itself
                        added[triI] = cellI;
                    }
                    else
                    {
                        added[triI] = meshMod.addCell
                        (
                            -1,     // masterPoint
                            -1,     // masterEdge
                            -1,     // masterFace
                            cellI,  // masterCell
                            mesh_.cellZones().whichZone(cellI)
                        );
                    }
                }
            }
        }
    }



    // Add triangle faces
    face triangle(3);

    forAll(mesh_.faces(), faceI)
    {
        label own = mesh_.faceOwner()[faceI];
        const labelList& addedOwn = faceOwnerCells_[faceI];
        const labelList& addedNei = faceNeighbourCells_[faceI];
        const face& f = mesh_.faces()[faceI];

        label patchI = -1;
        if (faceI >= mesh_.nInternalFaces())
        {
            patchI = mesh_.boundaryMesh().whichPatch(faceI);
        }

        label zoneI = mesh_.faceZones().whichZone(faceI);
        bool zoneFlip = false;
        if (zoneI != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            zoneFlip = fz.flipMap()[fz.whichFace(faceI)];
        }


        if (decomposeType == FACECENTRETETS)
        {
            forAll(f, fp)
            {
                // 1. Front triangle (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    triangle[0] = f[fp];
                    triangle[1] = f[f.fcIndex(fp)];
                    triangle[2] = faceToPoint_[faceI];

                    if (fp == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            triangle,
                            faceI,
                            addedOwn[fp],
                            addedNei[fp],
                            patchI,
                            zoneI,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            triangle,
                            addedOwn[fp],
                            addedNei[fp],
                            -1,                 //point
                            -1,                 //edge
                            faceI,              //face
                            patchI,
                            zoneI,
                            zoneFlip
                        );
                    }
                }


                // 2. Within owner cell - to cell centre
                {
                    label newOwn = addedOwn[f.rcIndex(fp)];
                    label newNei = addedOwn[fp];

                    triangle[0] = f[fp];
                    triangle[1] = cellToPoint_[own];
                    triangle[2] = faceToPoint_[faceI];

                    addFace
                    (
                        meshMod,
                        triangle,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchI
                        zoneI,
                        zoneFlip
                    );
                }
                // 2b. Within neighbour cell - to cell centre
                if (faceI < mesh_.nInternalFaces())
                {
                    label newOwn = addedNei[f.rcIndex(fp)];
                    label newNei = addedNei[fp];

                    triangle[0] = f[fp];
                    triangle[1] = faceToPoint_[faceI];
                    triangle[2] = cellToPoint_[mesh_.faceNeighbour()[faceI]];

                    addFace
                    (
                        meshMod,
                        triangle,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchI
                        zoneI,
                        zoneFlip
                    );
                }
            }
        }
        else
        {
            label fp0 = mesh_.tetBasePtIs()[faceI];
            label fp = f.fcIndex(fp0);

            for (label triI = 0; triI < f.size()-2; triI++)
            {
                label nextTri = triI+1;
                if (nextTri >= f.size()-2)
                {
                    nextTri -= f.size()-2;
                }
                label nextFp = f.fcIndex(fp);


                // Triangle triI consisiting of f[fp0], f[fp], f[nextFp]


                // 1. Front triangle (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    triangle[0] = f[fp0];
                    triangle[1] = f[fp];
                    triangle[2] = f[nextFp];

                    if (triI == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            triangle,
                            faceI,
                            addedOwn[triI],
                            addedNei[triI],
                            patchI,
                            zoneI,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            triangle,
                            addedOwn[triI],
                            addedNei[triI],
                            -1,                 //point
                            -1,                 //edge
                            faceI,              //face
                            patchI,
                            zoneI,
                            zoneFlip
                        );
                    }
                }


                // 2. Within owner cell - diagonal to cell centre
                if (triI < f.size()-3)
                {
                    label newOwn = addedOwn[triI];
                    label newNei = addedOwn[nextTri];

                    triangle[0] = f[fp0];
                    triangle[1] = f[nextFp];
                    triangle[2] = cellToPoint_[own];

                    addFace
                    (
                        meshMod,
                        triangle,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchI
                        zoneI,
                        zoneFlip
                    );

                    // 2b. Within neighbour cell - to cell centre
                    if (faceI < mesh_.nInternalFaces())
                    {
                        label newOwn = addedNei[triI];
                        label newNei = addedNei[nextTri];

                        triangle[0] = f[nextFp];
                        triangle[1] = f[fp0];
                        triangle[2] =
                            cellToPoint_[mesh_.faceNeighbour()[faceI]];

                        addFace
                        (
                            meshMod,
                            triangle,
                            newOwn,
                            newNei,
                            f[fp],      //point
                            -1,         //edge
                            -1,         //face
                            -1,         //patchI
                            zoneI,
                            zoneFlip
                        );
                    }
                }


                fp = nextFp;
            }
        }
    }



    // Add triangles for all edges.
    EdgeMap<label> edgeToFace;

    forAll(mesh_.cells(), cellI)
    {
        const cell& cFaces = mesh_.cells()[cellI];

        edgeToFace.clear();

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];

            label zoneI = mesh_.faceZones().whichZone(faceI);
            bool zoneFlip = false;
            if (zoneI != -1)
            {
                const faceZone& fz = mesh_.faceZones()[zoneI];
                zoneFlip = fz.flipMap()[fz.whichFace(faceI)];
            }

            const face& f = mesh_.faces()[faceI];
            //const labelList& fEdges = mesh_.faceEdges()[faceI];
            forAll(f, fp)
            {
                label p0 = f[fp];
                label p1 = f[f.fcIndex(fp)];
                const edge e(p0, p1);

                EdgeMap<label>::const_iterator edgeFnd = edgeToFace.find(e);
                if (edgeFnd == edgeToFace.end())
                {
                    edgeToFace.insert(e, faceI);
                }
                else
                {
                    // Found the other face on the edge.
                    label otherFaceI = edgeFnd();
                    const face& otherF = mesh_.faces()[otherFaceI];

                    // Found the other face on the edge. Note that since
                    // we are looping in the same order the tets added for
                    // otherFaceI will be before those of faceI

                    label otherFp = findIndex(otherF, p0);
                    if (otherF.nextLabel(otherFp) == p1)
                    {
                        // ok. otherFp is first vertex of edge.
                    }
                    else if (otherF.prevLabel(otherFp) == p1)
                    {
                        otherFp = otherF.rcIndex(otherFp);
                    }
                    else
                    {
                        FatalErrorIn("tetDecomposer::setRefinement(..)")
                            << "problem." << abort(FatalError);
                    }


                    // Triangle from edge to cell centre
                    if (mesh_.faceOwner()[faceI] == cellI)
                    {
                        triangle[0] = p0;
                        triangle[1] = p1;
                        triangle[2] = cellToPoint_[cellI];
                    }
                    else
                    {
                        triangle[0] = p1;
                        triangle[1] = p0;
                        triangle[2] = cellToPoint_[cellI];
                    }

                    // Determine tets on either side
                    label thisTet, otherTet;

                    if (decomposeType == FACECENTRETETS)
                    {
                        if (mesh_.faceOwner()[faceI] == cellI)
                        {
                            thisTet = faceOwnerCells_[faceI][fp];
                        }
                        else
                        {
                            thisTet = faceNeighbourCells_[faceI][fp];
                        }

                        if (mesh_.faceOwner()[otherFaceI] == cellI)
                        {
                            otherTet = faceOwnerCells_[otherFaceI][otherFp];
                        }
                        else
                        {
                            otherTet =
                                faceNeighbourCells_[otherFaceI][otherFp];
                        }
                    }
                    else
                    {
                        label thisTriI = triIndex(faceI, fp);
                        if (mesh_.faceOwner()[faceI] == cellI)
                        {
                            thisTet = faceOwnerCells_[faceI][thisTriI];
                        }
                        else
                        {
                            thisTet = faceNeighbourCells_[faceI][thisTriI];
                        }

                        label otherTriI = triIndex(otherFaceI, otherFp);
                        if (mesh_.faceOwner()[otherFaceI] == cellI)
                        {
                            otherTet = faceOwnerCells_[otherFaceI][otherTriI];
                        }
                        else
                        {
                            otherTet =
                                faceNeighbourCells_[otherFaceI][otherTriI];
                        }
                    }


                    addFace
                    (
                        meshMod,
                        triangle,
                        otherTet,
                        thisTet,
                        -1,         //masterPoint
                        -1,         //fEdges[fp], //masterEdge
                        faceI,      //masterFace
                        -1,         //patchI
                        zoneI,
                        zoneFlip
                    );
                }
            }
        }
    }
}


void Foam::tetDecomposer::updateMesh(const mapPolyMesh& map)
{
    inplaceRenumber(map.reversePointMap(), cellToPoint_);
    inplaceRenumber(map.reversePointMap(), faceToPoint_);

    forAll(faceOwnerCells_, faceI)
    {
        inplaceRenumber(map.reverseCellMap(), faceOwnerCells_[faceI]);
    }
    forAll(faceNeighbourCells_, faceI)
    {
        inplaceRenumber(map.reverseCellMap(), faceNeighbourCells_[faceI]);
    }
}


// ************************************************************************* //
