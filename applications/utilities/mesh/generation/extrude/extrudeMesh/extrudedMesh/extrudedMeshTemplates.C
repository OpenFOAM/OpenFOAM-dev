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

#include "extrudedMesh.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField
>
Foam::Xfer<Foam::pointField> Foam::extrudedMesh::extrudedPoints
(
    const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
    const extrudeModel& model
)
{
    const pointField& surfacePoints = extrudePatch.localPoints();
    const vectorField& surfaceNormals = extrudePatch.pointNormals();

    const label nLayers = model.nLayers();

    pointField ePoints((nLayers + 1)*surfacePoints.size());

    for (label layer=0; layer<=nLayers; layer++)
    {
        label offset = layer*surfacePoints.size();

        forAll(surfacePoints, i)
        {
            ePoints[offset + i] = model
            (
                surfacePoints[i],
                surfaceNormals[i],
                layer
            );
        }
    }

    // return points for transferring
    return xferMove(ePoints);
}


template<class Face, template<class> class FaceList, class PointField>
Foam::Xfer<Foam::faceList> Foam::extrudedMesh::extrudedFaces
(
    const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
    const extrudeModel& model
)
{
    const pointField& surfacePoints = extrudePatch.localPoints();
    const List<face>& surfaceFaces = extrudePatch.localFaces();
    const edgeList& surfaceEdges = extrudePatch.edges();
    const label nInternalEdges = extrudePatch.nInternalEdges();

    const label nLayers = model.nLayers();

    label nFaces =
        (nLayers + 1)*surfaceFaces.size() + nLayers*surfaceEdges.size();

    faceList eFaces(nFaces);

    labelList quad(4);
    label facei = 0;

    // Internal faces
    for (label layer=0; layer<nLayers; layer++)
    {
        label currentLayerOffset = layer*surfacePoints.size();
        label nextLayerOffset = currentLayerOffset + surfacePoints.size();

        // Vertical faces from layer to layer+1
        for (label edgeI=0; edgeI<nInternalEdges; edgeI++)
        {
            const edge& e = surfaceEdges[edgeI];
            const labelList& edgeFaces = extrudePatch.edgeFaces()[edgeI];

            face& f = eFaces[facei++];
            f.setSize(4);

            if
            (
                (edgeFaces[0] < edgeFaces[1])
             == sameOrder(surfaceFaces[edgeFaces[0]], e)
            )
            {
                f[0] = e[0] + currentLayerOffset;
                f[1] = e[1] + currentLayerOffset;
                f[2] = e[1] + nextLayerOffset;
                f[3] = e[0] + nextLayerOffset;
            }
            else
            {
                f[0] = e[1] + currentLayerOffset;
                f[1] = e[0] + currentLayerOffset;
                f[2] = e[0] + nextLayerOffset;
                f[3] = e[1] + nextLayerOffset;
            }
        }

        // Faces between layer and layer+1
        if (layer < nLayers-1)
        {
            forAll(surfaceFaces, i)
            {
                eFaces[facei++] =
                    face
                    (
                        surfaceFaces[i] //.reverseFace()
                      + nextLayerOffset
                    );
            }
        }
    }

    // External side faces
    for (label layer=0; layer<nLayers; layer++)
    {
        label currentLayerOffset = layer*surfacePoints.size();
        label nextLayerOffset = currentLayerOffset + surfacePoints.size();

        // Side faces across layer
        for (label edgeI=nInternalEdges; edgeI<surfaceEdges.size(); edgeI++)
        {
            const edge& e = surfaceEdges[edgeI];
            const labelList& edgeFaces = extrudePatch.edgeFaces()[edgeI];

            face& f = eFaces[facei++];
            f.setSize(4);

            if (sameOrder(surfaceFaces[edgeFaces[0]], e))
            {
                f[0] = e[0] + currentLayerOffset;
                f[1] = e[1] + currentLayerOffset;
                f[2] = e[1] + nextLayerOffset;
                f[3] = e[0] + nextLayerOffset;
            }
            else
            {
                f[0] = e[1] + currentLayerOffset;
                f[1] = e[0] + currentLayerOffset;
                f[2] = e[0] + nextLayerOffset;
                f[3] = e[1] + nextLayerOffset;
            }
        }
    }

    // Bottom faces
    forAll(surfaceFaces, i)
    {
        eFaces[facei++] = face(surfaceFaces[i]).reverseFace();
    }

    // Top faces
    forAll(surfaceFaces, i)
    {
        eFaces[facei++] =
            face
            (
                surfaceFaces[i]
              + nLayers*surfacePoints.size()
            );
    }

    // return points for transferring
    return xferMove(eFaces);
}


template<class Face, template<class> class FaceList, class PointField>
Foam::Xfer<Foam::cellList> Foam::extrudedMesh::extrudedCells
(
    const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
    const extrudeModel& model
)
{
    const List<face>& surfaceFaces = extrudePatch.localFaces();
    const edgeList& surfaceEdges = extrudePatch.edges();
    const label nInternalEdges = extrudePatch.nInternalEdges();

    const label nLayers = model.nLayers();

    cellList eCells(nLayers*surfaceFaces.size());

    // Size the cells
    forAll(surfaceFaces, i)
    {
        const face& f = surfaceFaces[i];

        for (label layer=0; layer<nLayers; layer++)
        {
            eCells[i + layer*surfaceFaces.size()].setSize(f.size() + 2);
        }
    }

    // Current face count per cell.
    labelList nCellFaces(eCells.size(), 0);


    label facei = 0;

    for (label layer=0; layer<nLayers; layer++)
    {
        // Side faces from layer to layer+1
        for (label i=0; i<nInternalEdges; i++)
        {
            // Get patch faces using edge
            const labelList& edgeFaces = extrudePatch.edgeFaces()[i];

            // Get cells on this layer
            label cell0 = layer*surfaceFaces.size() + edgeFaces[0];
            label cell1 = layer*surfaceFaces.size() + edgeFaces[1];

            eCells[cell0][nCellFaces[cell0]++] = facei;
            eCells[cell1][nCellFaces[cell1]++] = facei;

            facei++;
        }

        // Faces between layer and layer+1
        if (layer < nLayers-1)
        {
            forAll(surfaceFaces, i)
            {
                label cell0 = layer*surfaceFaces.size() + i;
                label cell1 = (layer+1)*surfaceFaces.size() + i;

                eCells[cell0][nCellFaces[cell0]++] = facei;
                eCells[cell1][nCellFaces[cell1]++] = facei;

                facei++;
            }
        }
    }

    // External side faces
    for (label layer=0; layer<nLayers; layer++)
    {
        // Side faces across layer
        for (label i=nInternalEdges; i<surfaceEdges.size(); i++)
        {
            // Get patch faces using edge
            const labelList& edgeFaces = extrudePatch.edgeFaces()[i];

            // Get cells on this layer
            label cell0 = layer*surfaceFaces.size() + edgeFaces[0];

            eCells[cell0][nCellFaces[cell0]++] = facei;

            facei++;
        }
    }

    // Top faces
    forAll(surfaceFaces, i)
    {
        eCells[i][nCellFaces[i]++] = facei;

        facei++;
    }

    // Bottom faces
    forAll(surfaceFaces, i)
    {
        label cell0 = (nLayers-1)*surfaceFaces.size() + i;

        eCells[cell0][nCellFaces[cell0]++] = facei;

        facei++;
    }

    // return points for transferring
    return xferMove(eCells);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField
>
Foam::extrudedMesh::extrudedMesh
(
    const IOobject& io,
    const PrimitivePatch<Face, FaceList, PointField>& extrudePatch,
    const extrudeModel& model
)
:
    polyMesh
    (
        io,
        extrudedPoints(extrudePatch, model),
        extrudedFaces(extrudePatch, model),
        extrudedCells(extrudePatch, model)
    ),
    model_(model)
{
    List<polyPatch*> patches(3);

    label facei = nInternalFaces();

    label sz =
        model_.nLayers()
       *(extrudePatch.nEdges() - extrudePatch.nInternalEdges());

    patches[0] = new wallPolyPatch
    (
        "sides",
        sz,
        facei,
        0,
        boundaryMesh(),
        wallPolyPatch::typeName
    );

    facei += sz;

    patches[1] = new polyPatch
    (
        "originalPatch",
        extrudePatch.size(),
        facei,
        1,
        boundaryMesh(),
        polyPatch::typeName
    );

    facei += extrudePatch.size();

    patches[2] = new polyPatch
    (
        "otherSide",
        extrudePatch.size(),
        facei,
        2,
        boundaryMesh(),
        polyPatch::typeName
    );

    addPatches(patches);
}


// ************************************************************************* //
