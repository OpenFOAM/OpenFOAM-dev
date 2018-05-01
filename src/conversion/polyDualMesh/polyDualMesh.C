/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

InClass
    polyDualMesh

\*---------------------------------------------------------------------------*/

#include "polyDualMesh.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Time.H"
#include "SortableList.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyDualMesh, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Determine order for faces:
// - upper-triangular order for internal faces
// - external faces after internal faces (were already so)
Foam::labelList Foam::polyDualMesh::getFaceOrder
(
    const labelList& faceOwner,
    const labelList& faceNeighbour,
    const cellList& cells,
    label& nInternalFaces
)
{
    labelList oldToNew(faceOwner.size(), -1);

    // First unassigned face
    label newFacei = 0;

    forAll(cells, celli)
    {
        const labelList& cFaces = cells[celli];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            label nbrCelli = faceNeighbour[facei];

            if (nbrCelli != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCelli == celli)
                {
                    nbrCelli = faceOwner[facei];
                }

                if (celli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNew[cFaces[nbr.indices()[i]]] = newFacei++;
            }
        }
    }

    nInternalFaces = newFacei;

    Pout<< "nInternalFaces:" << nInternalFaces << endl;
    Pout<< "nFaces:" << faceOwner.size() << endl;
    Pout<< "nCells:" << cells.size() << endl;

    // Leave patch faces intact.
    for (label facei = newFacei; facei < faceOwner.size(); facei++)
    {
        oldToNew[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNew, facei)
    {
        if (oldToNew[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position"
                << " for face " << facei
                << abort(FatalError);
        }
    }

    return oldToNew;
}


// Get the two edges on facei using pointi. Returns them such that the order
// - otherVertex of e0
// - pointi
// - otherVertex(pointi) of e1
// is in face order
void Foam::polyDualMesh::getPointEdges
(
    const primitivePatch& patch,
    const label facei,
    const label pointi,
    label& e0,
    label& e1
)
{
    const labelList& fEdges = patch.faceEdges()[facei];
    const face& f = patch.localFaces()[facei];

    e0 = -1;
    e1 = -1;

    forAll(fEdges, i)
    {
        label edgeI = fEdges[i];

        const edge& e = patch.edges()[edgeI];

        if (e[0] == pointi)
        {
            // One of the edges using pointi. Check which one.
            label index = findIndex(f, pointi);

            if (f.nextLabel(index) == e[1])
            {
                e1 = edgeI;
            }
            else
            {
                e0 = edgeI;
            }

            if (e0 != -1 && e1 != -1)
            {
                return;
            }
        }
        else if (e[1] == pointi)
        {
            // One of the edges using pointi. Check which one.
            label index = findIndex(f, pointi);

            if (f.nextLabel(index) == e[0])
            {
                e1 = edgeI;
            }
            else
            {
                e0 = edgeI;
            }

            if (e0 != -1 && e1 != -1)
            {
                return;
            }
        }
    }

    FatalErrorInFunction
        << " vertices:" << patch.localFaces()[facei]
        << " that uses point:" << pointi
        << abort(FatalError);
}


// Collect the face on an external point of the patch.
Foam::labelList Foam::polyDualMesh::collectPatchSideFace
(
    const polyPatch& patch,
    const label patchToDualOffset,
    const labelList& edgeToDualPoint,
    const labelList& pointToDualPoint,
    const label pointi,

    label& edgeI
)
{
    // Construct face by walking around e[eI] starting from
    // patchEdgeI

    label meshPointi = patch.meshPoints()[pointi];
    const labelList& pFaces = patch.pointFaces()[pointi];

    DynamicList<label> dualFace;

    if (pointToDualPoint[meshPointi] >= 0)
    {
        // Number of pFaces + 2 boundary edge + feature point
        dualFace.setCapacity(pFaces.size()+2+1);
        // Store dualVertex for feature edge
        dualFace.append(pointToDualPoint[meshPointi]);
    }
    else
    {
        dualFace.setCapacity(pFaces.size()+2);
    }

    // Store dual vertex for starting edge.
    if (edgeToDualPoint[patch.meshEdges()[edgeI]] < 0)
    {
        FatalErrorInFunction
            << abort(FatalError);
    }

    dualFace.append(edgeToDualPoint[patch.meshEdges()[edgeI]]);

    label facei = patch.edgeFaces()[edgeI][0];

    // Check order of vertices of edgeI in face to see if we need to reverse.
    bool reverseFace;

    label e0, e1;
    getPointEdges(patch, facei, pointi, e0, e1);

    if (e0 == edgeI)
    {
        reverseFace = true;
    }
    else
    {
        reverseFace = false;
    }

    while (true)
    {
        // Store dual vertex for facei.
        dualFace.append(facei + patchToDualOffset);

        // Cross face to other edge on pointi
        label e0, e1;
        getPointEdges(patch, facei, pointi, e0, e1);

        if (e0 == edgeI)
        {
            edgeI = e1;
        }
        else
        {
            edgeI = e0;
        }

        if (edgeToDualPoint[patch.meshEdges()[edgeI]] >= 0)
        {
            // Feature edge. Insert dual vertex for edge.
            dualFace.append(edgeToDualPoint[patch.meshEdges()[edgeI]]);
        }

        const labelList& eFaces = patch.edgeFaces()[edgeI];

        if (eFaces.size() == 1)
        {
            // Reached other edge of patch
            break;
        }

        // Cross edge to other face.
        if (eFaces[0] == facei)
        {
            facei = eFaces[1];
        }
        else
        {
            facei = eFaces[0];
        }
    }

    dualFace.shrink();

    if (reverseFace)
    {
        reverse(dualFace);
    }

    return dualFace;
}


// Collect face around pointi which is not on the outside of the patch.
// Returns the vertices of the face and the indices in these vertices of
// any points which are on feature edges.
void Foam::polyDualMesh::collectPatchInternalFace
(
    const polyPatch& patch,
    const label patchToDualOffset,
    const labelList& edgeToDualPoint,

    const label pointi,
    const label startEdgeI,

    labelList& dualFace2,
    labelList& featEdgeIndices2
)
{
    // Construct face by walking around pointi starting from startEdgeI
    const labelList& meshEdges = patch.meshEdges();
    const labelList& pFaces = patch.pointFaces()[pointi];

    // Vertices of dualFace
    DynamicList<label> dualFace(pFaces.size());
    // Indices in dualFace of vertices added because of feature edge
    DynamicList<label> featEdgeIndices(dualFace.size());


    label edgeI = startEdgeI;
    label facei = patch.edgeFaces()[edgeI][0];

    // Check order of vertices of edgeI in face to see if we need to reverse.
    bool reverseFace;

    label e0, e1;
    getPointEdges(patch, facei, pointi, e0, e1);

    if (e0 == edgeI)
    {
        reverseFace = true;
    }
    else
    {
        reverseFace = false;
    }

    while (true)
    {
        // Insert dual vertex for face
        dualFace.append(facei + patchToDualOffset);

        // Cross face to other edge on pointi
        label e0, e1;
        getPointEdges(patch, facei, pointi, e0, e1);

        if (e0 == edgeI)
        {
            edgeI = e1;
        }
        else
        {
            edgeI = e0;
        }

        if (edgeToDualPoint[meshEdges[edgeI]] >= 0)
        {
            // Feature edge. Insert dual vertex for edge.
            dualFace.append(edgeToDualPoint[meshEdges[edgeI]]);
            featEdgeIndices.append(dualFace.size()-1);
        }

        if (edgeI == startEdgeI)
        {
            break;
        }

        // Cross edge to other face.
        const labelList& eFaces = patch.edgeFaces()[edgeI];

        if (eFaces[0] == facei)
        {
            facei = eFaces[1];
        }
        else
        {
            facei = eFaces[0];
        }
    }

    dualFace2.transfer(dualFace);

    featEdgeIndices2.transfer(featEdgeIndices);

    if (reverseFace)
    {
        reverse(dualFace2);

        // Correct featEdgeIndices for change in dualFace2
        forAll(featEdgeIndices2, i)
        {
            featEdgeIndices2[i] = dualFace2.size() -1 - featEdgeIndices2[i];
        }
        // Reverse indices (might not be necessary but do anyway)
        reverse(featEdgeIndices2);
    }
}


void Foam::polyDualMesh::splitFace
(
    const polyPatch& patch,
    const labelList& pointToDualPoint,

    const label pointi,
    const labelList& dualFace,
    const labelList& featEdgeIndices,

    DynamicList<face>& dualFaces,
    DynamicList<label>& dualOwner,
    DynamicList<label>& dualNeighbour,
    DynamicList<label>& dualRegion
)
{

    // Split face because of feature edges/points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label meshPointi = patch.meshPoints()[pointi];

    if (pointToDualPoint[meshPointi] >= 0)
    {
        // Feature point. Do face-centre decomposition.

        if (featEdgeIndices.size() < 2)
        {
            // Feature point but no feature edges. Not handled at the moment
            dualFaces.append(face(dualFace));
            dualOwner.append(meshPointi);
            dualNeighbour.append(-1);
            dualRegion.append(patch.index());
        }
        else
        {
            // Do 'face-centre' decomposition. Start from first feature
            // edge create face up until next feature edge.

            forAll(featEdgeIndices, i)
            {
                label startFp = featEdgeIndices[i];

                label endFp = featEdgeIndices[(i+1) % featEdgeIndices.size()];

                // Collect face from startFp to endFp
                label sz = 0;

                if (endFp > startFp)
                {
                    sz = endFp - startFp + 2;
                }
                else
                {
                    sz = endFp + dualFace.size() - startFp + 2;
                }
                face subFace(sz);

                // feature point becomes face centre.
                subFace[0] = pointToDualPoint[patch.meshPoints()[pointi]];

                // Copy from startFp up to endFp.
                for (label subFp = 1; subFp < subFace.size(); subFp++)
                {
                    subFace[subFp] = dualFace[startFp];

                    startFp = (startFp+1) % dualFace.size();
                }

                dualFaces.append(face(subFace));
                dualOwner.append(meshPointi);
                dualNeighbour.append(-1);
                dualRegion.append(patch.index());
            }
        }
    }
    else
    {
        // No feature point. Check if feature edges.
        if (featEdgeIndices.size() < 2)
        {
            // Not enough feature edges. No split.
            dualFaces.append(face(dualFace));
            dualOwner.append(meshPointi);
            dualNeighbour.append(-1);
            dualRegion.append(patch.index());
        }
        else
        {
            // Multiple feature edges but no feature point.
            // Split face along feature edges. Gives weird result if
            // number of feature edges > 2.

            // Storage for new face
            DynamicList<label> subFace(dualFace.size());

            forAll(featEdgeIndices, featI)
            {
                label startFp = featEdgeIndices[featI];
                label endFp = featEdgeIndices[featEdgeIndices.fcIndex(featI)];

                label fp = startFp;

                while (true)
                {
                    label vertI = dualFace[fp];

                    subFace.append(vertI);

                    if (fp == endFp)
                    {
                        break;
                    }

                    fp = dualFace.fcIndex(fp);
                }

                if (subFace.size() > 2)
                {
                    // Enough vertices to create a face from.
                    subFace.shrink();

                    dualFaces.append(face(subFace));
                    dualOwner.append(meshPointi);
                    dualNeighbour.append(-1);
                    dualRegion.append(patch.index());

                    subFace.clear();
                }
            }
            // Check final face.
            if (subFace.size() > 2)
            {
                // Enough vertices to create a face from.
                subFace.shrink();

                dualFaces.append(face(subFace));
                dualOwner.append(meshPointi);
                dualNeighbour.append(-1);
                dualRegion.append(patch.index());

                subFace.clear();
            }
        }
    }
}


// Create boundary face for every point in patch
void Foam::polyDualMesh::dualPatch
(
    const polyPatch& patch,
    const label patchToDualOffset,
    const labelList& edgeToDualPoint,
    const labelList& pointToDualPoint,

    const pointField& dualPoints,

    DynamicList<face>& dualFaces,
    DynamicList<label>& dualOwner,
    DynamicList<label>& dualNeighbour,
    DynamicList<label>& dualRegion
)
{
    const labelList& meshEdges = patch.meshEdges();

    // Whether edge has been done.
    // 0 : not
    // 1 : done e.start()
    // 2 : done e.end()
    // 3 : done both
    labelList doneEdgeSide(meshEdges.size(), 0);

    boolList donePoint(patch.nPoints(), false);


    // Do points on edge of patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(doneEdgeSide, patchEdgeI)
    {
        const labelList& eFaces = patch.edgeFaces()[patchEdgeI];

        if (eFaces.size() == 1)
        {
            const edge& e = patch.edges()[patchEdgeI];

            forAll(e, eI)
            {
                label bitMask = 1<<eI;

                if ((doneEdgeSide[patchEdgeI] & bitMask) == 0)
                {
                    // Construct face by walking around e[eI] starting from
                    // patchEdgeI

                    label pointi = e[eI];

                    label edgeI = patchEdgeI;
                    labelList dualFace
                    (
                        collectPatchSideFace
                        (
                            patch,
                            patchToDualOffset,
                            edgeToDualPoint,
                            pointToDualPoint,

                            pointi,
                            edgeI
                        )
                    );

                    // Now edgeI is end edge. Mark as visited
                    if (patch.edges()[edgeI][0] == pointi)
                    {
                        doneEdgeSide[edgeI] |= 1;
                    }
                    else
                    {
                        doneEdgeSide[edgeI] |= 2;
                    }

                    dualFaces.append(face(dualFace));
                    dualOwner.append(patch.meshPoints()[pointi]);
                    dualNeighbour.append(-1);
                    dualRegion.append(patch.index());

                    doneEdgeSide[patchEdgeI] |= bitMask;
                    donePoint[pointi] = true;
                }
            }
        }
    }



    // Do patch-internal points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(donePoint, pointi)
    {
        if (!donePoint[pointi])
        {
            labelList dualFace, featEdgeIndices;

            collectPatchInternalFace
            (
                patch,
                patchToDualOffset,
                edgeToDualPoint,
                pointi,
                patch.pointEdges()[pointi][0],  // Arbitrary starting edge

                dualFace,
                featEdgeIndices
            );

            //- Either keep in one piece or split face according to feature.

            //// Keep face in one piece.
            // dualFaces.append(face(dualFace));
            // dualOwner.append(patch.meshPoints()[pointi]);
            // dualNeighbour.append(-1);
            // dualRegion.append(patch.index());

            splitFace
            (
                patch,
                pointToDualPoint,
                pointi,
                dualFace,
                featEdgeIndices,

                dualFaces,
                dualOwner,
                dualNeighbour,
                dualRegion
            );

            donePoint[pointi] = true;
        }
    }
}


void Foam::polyDualMesh::calcDual
(
    const polyMesh& mesh,
    const labelList& featureEdges,
    const labelList& featurePoints
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const label nIntFaces = mesh.nInternalFaces();


    // Get patch for all of outside
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces() - nIntFaces,
            nIntFaces
        ),
        mesh.points()
    );

    // Correspondence from patch edge to mesh edge.
    labelList meshEdges
    (
        allBoundary.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    {
        pointSet nonManifoldPoints
        (
            mesh,
            "nonManifoldPoints",
            meshEdges.size() / 100
        );

        allBoundary.checkPointManifold(true, &nonManifoldPoints);

        if (nonManifoldPoints.size())
        {
            nonManifoldPoints.write();

            FatalErrorInFunction
                << "There are " << nonManifoldPoints.size() << " points where"
                << " the outside of the mesh is non-manifold." << nl
                << "Such a mesh cannot be converted to a dual." << nl
                << "Writing points at non-manifold sites to pointSet "
                << nonManifoldPoints.name()
                << exit(FatalError);
        }
    }


    // Assign points
    // ~~~~~~~~~~~~~

    // mesh label   dualMesh vertex
    // ----------   ---------------
    // celli        celli
    // facei        nCells+facei-nIntFaces
    // featEdgeI    nCells+nFaces-nIntFaces+featEdgeI
    // featPointi   nCells+nFaces-nIntFaces+nFeatEdges+featPointi

    pointField dualPoints
    (
        mesh.nCells()                           // cell centres
      + mesh.nFaces() - nIntFaces               // boundary face centres
      + featureEdges.size()                     // additional boundary edges
      + featurePoints.size()                    // additional boundary points
    );

    label dualPointi = 0;


    // Cell centres.
    const pointField& cellCentres = mesh.cellCentres();

    cellPoint_.setSize(cellCentres.size());

    forAll(cellCentres, celli)
    {
        cellPoint_[celli] = dualPointi;
        dualPoints[dualPointi++] = cellCentres[celli];
    }

    // Boundary faces centres
    const pointField& faceCentres = mesh.faceCentres();

    boundaryFacePoint_.setSize(mesh.nFaces() - nIntFaces);

    for (label facei = nIntFaces; facei < mesh.nFaces(); facei++)
    {
        boundaryFacePoint_[facei - nIntFaces] = dualPointi;
        dualPoints[dualPointi++] = faceCentres[facei];
    }

    // Edge status:
    //  >0 : dual point label (edge is feature edge)
    //  -1 : is boundary edge.
    //  -2 : is internal edge.
    labelList edgeToDualPoint(mesh.nEdges(), -2);

    forAll(meshEdges, patchEdgeI)
    {
        label edgeI = meshEdges[patchEdgeI];

        edgeToDualPoint[edgeI] = -1;
    }

    forAll(featureEdges, i)
    {
        label edgeI = featureEdges[i];

        const edge& e = mesh.edges()[edgeI];

        edgeToDualPoint[edgeI] = dualPointi;
        dualPoints[dualPointi++] = e.centre(mesh.points());
    }



    // Point status:
    //  >0 : dual point label (point is feature point)
    //  -1 : is point on edge of patch
    //  -2 : is point on patch (but not on edge)
    //  -3 : is internal point.
    labelList pointToDualPoint(mesh.nPoints(), -3);

    forAll(patches, patchi)
    {
        const labelList& meshPoints = patches[patchi].meshPoints();

        forAll(meshPoints, i)
        {
            pointToDualPoint[meshPoints[i]] = -2;
        }

        const labelListList& loops = patches[patchi].edgeLoops();

        forAll(loops, i)
        {
            const labelList& loop = loops[i];

            forAll(loop, j)
            {
                 pointToDualPoint[meshPoints[loop[j]]] = -1;
            }
        }
    }

    forAll(featurePoints, i)
    {
        label pointi = featurePoints[i];

        pointToDualPoint[pointi] = dualPointi;
        dualPoints[dualPointi++] = mesh.points()[pointi];
    }


    // Storage for new faces.
    // Dynamic sized since we don't know size.

    DynamicList<face> dynDualFaces(mesh.nEdges());
    DynamicList<label> dynDualOwner(mesh.nEdges());
    DynamicList<label> dynDualNeighbour(mesh.nEdges());
    DynamicList<label> dynDualRegion(mesh.nEdges());


    // Generate faces from edges on the boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(meshEdges, patchEdgeI)
    {
        label edgeI = meshEdges[patchEdgeI];

        const edge& e = mesh.edges()[edgeI];

        label owner = -1;
        label neighbour = -1;

        if (e[0] < e[1])
        {
            owner = e[0];
            neighbour = e[1];
        }
        else
        {
            owner = e[1];
            neighbour = e[0];
        }

        // Find the boundary faces using the edge.
        const labelList& patchFaces = allBoundary.edgeFaces()[patchEdgeI];

        if (patchFaces.size() != 2)
        {
            FatalErrorInFunction
                << "Cannot handle edges with " << patchFaces.size()
                << " connected boundary faces."
                << abort(FatalError);
        }

        label face0 = patchFaces[0] + nIntFaces;
        const face& f0 = mesh.faces()[face0];

        label face1 = patchFaces[1] + nIntFaces;

        // We want to start walking from patchFaces[0] or patchFaces[1],
        // depending on which one uses owner,neighbour in the right order.

        label startFacei = -1;
        label endFacei = -1;

        label index = findIndex(f0, neighbour);

        if (f0.nextLabel(index) == owner)
        {
            startFacei = face0;
            endFacei = face1;
        }
        else
        {
            startFacei = face1;
            endFacei = face0;
        }

        // Now walk from startFacei to cell to face to cell etc. until endFacei

        DynamicList<label> dualFace;

        if (edgeToDualPoint[edgeI] >= 0)
        {
            // Number of cells + 2 boundary faces + feature edge point
            dualFace.setCapacity(mesh.edgeCells()[edgeI].size()+2+1);
            // Store dualVertex for feature edge
            dualFace.append(edgeToDualPoint[edgeI]);
        }
        else
        {
            dualFace.setCapacity(mesh.edgeCells()[edgeI].size()+2);
        }

        // Store dual vertex for starting face.
        dualFace.append(mesh.nCells() + startFacei - nIntFaces);

        label celli = mesh.faceOwner()[startFacei];
        label facei = startFacei;

        while (true)
        {
            // Store dual vertex from celli.
            dualFace.append(celli);

            // Cross cell to other face on edge.
            label f0, f1;
            meshTools::getEdgeFaces(mesh, celli, edgeI, f0, f1);

            if (f0 == facei)
            {
                facei = f1;
            }
            else
            {
                facei = f0;
            }

            // Cross face to other cell.
            if (facei == endFacei)
            {
                break;
            }

            if (mesh.faceOwner()[facei] == celli)
            {
                celli = mesh.faceNeighbour()[facei];
            }
            else
            {
                celli = mesh.faceOwner()[facei];
            }
        }

        // Store dual vertex for endFace.
        dualFace.append(mesh.nCells() + endFacei - nIntFaces);

        dynDualFaces.append(face(dualFace.shrink()));
        dynDualOwner.append(owner);
        dynDualNeighbour.append(neighbour);
        dynDualRegion.append(-1);

        {
            // Check orientation.
            const face& f = dynDualFaces.last();
            const vector a = f.area(dualPoints);
            if (((mesh.points()[owner] - dualPoints[f[0]]) & a) > 0)
            {
                WarningInFunction
                    << " on boundary edge:" << edgeI
                    << mesh.points()[mesh.edges()[edgeI][0]]
                    << mesh.points()[mesh.edges()[edgeI][1]]
                    << endl;
            }
        }
    }


    // Generate faces from internal edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(edgeToDualPoint, edgeI)
    {
        if (edgeToDualPoint[edgeI] == -2)
        {
            // Internal edge.

            // Find dual owner, neighbour

            const edge& e = mesh.edges()[edgeI];

            label owner = -1;
            label neighbour = -1;

            if (e[0] < e[1])
            {
                owner = e[0];
                neighbour = e[1];
            }
            else
            {
                owner = e[1];
                neighbour = e[0];
            }

            // Get a starting cell
            const labelList& eCells = mesh.edgeCells()[edgeI];

            label celli = eCells[0];

            // Get the two faces on the cell and edge.
            label face0, face1;
            meshTools::getEdgeFaces(mesh, celli, edgeI, face0, face1);

            // Find the starting face by looking at the order in which
            // the face uses the owner, neighbour
            const face& f0 = mesh.faces()[face0];

            label index = findIndex(f0, neighbour);

            bool f0OrderOk = (f0.nextLabel(index) == owner);

            label startFacei = -1;

            if (f0OrderOk == (mesh.faceOwner()[face0] == celli))
            {
                startFacei = face0;
            }
            else
            {
                startFacei = face1;
            }

            // Walk face-cell-face until starting face reached.
            DynamicList<label> dualFace(mesh.edgeCells()[edgeI].size());

            label facei = startFacei;

            while (true)
            {
                // Store dual vertex from celli.
                dualFace.append(celli);

                // Cross cell to other face on edge.
                label f0, f1;
                meshTools::getEdgeFaces(mesh, celli, edgeI, f0, f1);

                if (f0 == facei)
                {
                    facei = f1;
                }
                else
                {
                    facei = f0;
                }

                // Cross face to other cell.
                if (facei == startFacei)
                {
                    break;
                }

                if (mesh.faceOwner()[facei] == celli)
                {
                    celli = mesh.faceNeighbour()[facei];
                }
                else
                {
                    celli = mesh.faceOwner()[facei];
                }
            }

            dynDualFaces.append(face(dualFace.shrink()));
            dynDualOwner.append(owner);
            dynDualNeighbour.append(neighbour);
            dynDualRegion.append(-1);

            {
                // Check orientation.
                const face& f = dynDualFaces.last();
                const vector a = f.area(dualPoints);
                if (((mesh.points()[owner] - dualPoints[f[0]]) & a) > 0)
                {
                    WarningInFunction
                        << " on internal edge:" << edgeI
                        << mesh.points()[mesh.edges()[edgeI][0]]
                        << mesh.points()[mesh.edges()[edgeI][1]]
                        << endl;
                }
            }
        }
    }

    // Dump faces.
    if (debug)
    {
        dynDualFaces.shrink();
        dynDualOwner.shrink();
        dynDualNeighbour.shrink();
        dynDualRegion.shrink();

        OFstream str("dualInternalFaces.obj");
        Pout<< "polyDualMesh::calcDual : dumping internal faces to "
            << str.name() << endl;

        forAll(dualPoints, dualPointi)
        {
            meshTools::writeOBJ(str, dualPoints[dualPointi]);
        }
        forAll(dynDualFaces, dualFacei)
        {
            const face& f = dynDualFaces[dualFacei];

            str<< 'f';
            forAll(f, fp)
            {
                str<< ' ' << f[fp]+1;
            }
            str<< nl;
        }
    }

    const label nInternalFaces = dynDualFaces.size();

    // Outside faces
    // ~~~~~~~~~~~~~

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        dualPatch
        (
            pp,
            mesh.nCells() + pp.start() - nIntFaces,
            edgeToDualPoint,
            pointToDualPoint,

            dualPoints,

            dynDualFaces,
            dynDualOwner,
            dynDualNeighbour,
            dynDualRegion
        );
    }


    // Transfer face info to straight lists
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    faceList dualFaces(dynDualFaces.shrink(), true);
    dynDualFaces.clear();

    labelList dualOwner(dynDualOwner.shrink(), true);
    dynDualOwner.clear();

    labelList dualNeighbour(dynDualNeighbour.shrink(), true);
    dynDualNeighbour.clear();

    labelList dualRegion(dynDualRegion.shrink(), true);
    dynDualRegion.clear();



    // Dump faces.
    if (debug)
    {
        OFstream str("dualFaces.obj");
        Pout<< "polyDualMesh::calcDual : dumping all faces to "
            << str.name() << endl;

        forAll(dualPoints, dualPointi)
        {
            meshTools::writeOBJ(str, dualPoints[dualPointi]);
        }
        forAll(dualFaces, dualFacei)
        {
            const face& f = dualFaces[dualFacei];

            str<< 'f';
            forAll(f, fp)
            {
                str<< ' ' << f[fp]+1;
            }
            str<< nl;
        }
    }

    // Create cells.
    cellList dualCells(mesh.nPoints());

    forAll(dualCells, celli)
    {
        dualCells[celli].setSize(0);
    }

    forAll(dualOwner, facei)
    {
        label celli = dualOwner[facei];

        labelList& cFaces = dualCells[celli];

        label sz = cFaces.size();
        cFaces.setSize(sz+1);
        cFaces[sz] = facei;
    }
    forAll(dualNeighbour, facei)
    {
        label celli = dualNeighbour[facei];

        if (celli != -1)
        {
            labelList& cFaces = dualCells[celli];

            label sz = cFaces.size();
            cFaces.setSize(sz+1);
            cFaces[sz] = facei;
        }
    }


    // Do upper-triangular face ordering. Determines face reordering map and
    // number of internal faces.
    label dummy;

    labelList oldToNew
    (
        getFaceOrder
        (
            dualOwner,
            dualNeighbour,
            dualCells,
            dummy
        )
    );

    // Reorder faces.
    inplaceReorder(oldToNew, dualFaces);
    inplaceReorder(oldToNew, dualOwner);
    inplaceReorder(oldToNew, dualNeighbour);
    inplaceReorder(oldToNew, dualRegion);
    forAll(dualCells, celli)
    {
        inplaceRenumber(oldToNew, dualCells[celli]);
    }


    // Create patches
    labelList patchSizes(patches.size(), 0);

    forAll(dualRegion, facei)
    {
        if (dualRegion[facei] >= 0)
        {
            patchSizes[dualRegion[facei]]++;
        }
    }

    labelList patchStarts(patches.size(), 0);

    label facei = nInternalFaces;

    forAll(patches, patchi)
    {
        patchStarts[patchi] = facei;

        facei += patchSizes[patchi];
    }


    Pout<< "nFaces:" << dualFaces.size()
        << " patchSizes:" << patchSizes
        << " patchStarts:" << patchStarts
        << endl;


    // Add patches. First add zero sized (since mesh still 0 size)
    List<polyPatch*> dualPatches(patches.size());

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        dualPatches[patchi] = pp.clone
        (
            boundaryMesh(),
            patchi,
            0, // patchSizes[patchi],
            0  // patchStarts[patchi]
        ).ptr();
    }
    addPatches(dualPatches);

    // Assign to mesh.
    resetPrimitives
    (
        xferMove(dualPoints),
        xferMove(dualFaces),
        xferMove(dualOwner),
        xferMove(dualNeighbour),
        patchSizes,
        patchStarts
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyDualMesh::polyDualMesh(const IOobject& io)
:
    polyMesh(io),
    cellPoint_
    (
        IOobject
        (
            "cellPoint",
            time().findInstance(meshDir(), "cellPoint"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundaryFacePoint_
    (
        IOobject
        (
            "boundaryFacePoint",
            time().findInstance(meshDir(), "boundaryFacePoint"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{}


// Construct from polyMesh
Foam::polyDualMesh::polyDualMesh
(
    const polyMesh& mesh,
    const labelList& featureEdges,
    const labelList& featurePoints
)
:
    polyMesh
    (
        mesh,
        xferCopy(pointField()),// to prevent any warnings "points not allocated"
        xferCopy(faceList()),  // to prevent any warnings "faces  not allocated"
        xferCopy(cellList())
    ),
    cellPoint_
    (
        IOobject
        (
            "cellPoint",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nCells(), -1)
    ),
    boundaryFacePoint_
    (
        IOobject
        (
            "boundaryFacePoint",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nFaces() - mesh.nInternalFaces())
    )
{
    calcDual(mesh, featureEdges, featurePoints);
}


// Construct from polyMesh and feature angle
Foam::polyDualMesh::polyDualMesh
(
    const polyMesh& mesh,
    const scalar featureCos
)
:
    polyMesh
    (
        mesh,
        xferCopy(pointField()),// to prevent any warnings "points not allocated"
        xferCopy(faceList()),  // to prevent any warnings "faces  not allocated"
        xferCopy(cellList())
    ),
    cellPoint_
    (
        IOobject
        (
            "cellPoint",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nCells(), -1)
    ),
    boundaryFacePoint_
    (
        IOobject
        (
            "boundaryFacePoint",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nFaces() - mesh.nInternalFaces(), -1)
    )
{
    labelList featureEdges, featurePoints;

    calcFeatures(mesh, featureCos, featureEdges, featurePoints);
    calcDual(mesh, featureEdges, featurePoints);
}


void Foam::polyDualMesh::calcFeatures
(
    const polyMesh& mesh,
    const scalar featureCos,
    labelList& featureEdges,
    labelList& featurePoints
)
{
    // Create big primitivePatch for all outside.
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces() - mesh.nInternalFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // For ease of use store patch number per face in allBoundary.
    labelList allRegion(allBoundary.size());

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, i)
        {
            allRegion[i + pp.start() - mesh.nInternalFaces()] = patchi;
        }
    }


    // Calculate patch/feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelListList& edgeFaces = allBoundary.edgeFaces();
    const vectorField& faceNormals = allBoundary.faceNormals();
    const labelList& meshPoints = allBoundary.meshPoints();

    boolList isFeatureEdge(edgeFaces.size(), false);

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() != 2)
        {
            // Non-manifold. Problem?
            const edge& e = allBoundary.edges()[edgeI];

            WarningInFunction
                << meshPoints[e[0]] << ' ' << meshPoints[e[1]]
                << "  coords:" << mesh.points()[meshPoints[e[0]]]
                << mesh.points()[meshPoints[e[1]]]
                << " has more than 2 faces connected to it:"
                << eFaces.size() << endl;

            isFeatureEdge[edgeI] = true;
        }
        else if (allRegion[eFaces[0]] != allRegion[eFaces[1]])
        {
            isFeatureEdge[edgeI] = true;
        }
        else if
        (
            (faceNormals[eFaces[0]] & faceNormals[eFaces[1]])
          < featureCos
        )
        {
            isFeatureEdge[edgeI] = true;
        }
    }


    // Calculate feature points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    const labelListList& pointEdges = allBoundary.pointEdges();

    DynamicList<label> allFeaturePoints(pointEdges.size());

    forAll(pointEdges, pointi)
    {
        const labelList& pEdges = pointEdges[pointi];

        label nFeatEdges = 0;

        forAll(pEdges, i)
        {
            if (isFeatureEdge[pEdges[i]])
            {
                nFeatEdges++;
            }
        }
        if (nFeatEdges > 2)
        {
            // Store in mesh vertex label
            allFeaturePoints.append(allBoundary.meshPoints()[pointi]);
        }
    }
    featurePoints.transfer(allFeaturePoints);

    if (debug)
    {
        OFstream str("featurePoints.obj");

        Pout<< "polyDualMesh::calcFeatures : dumping feature points to "
            << str.name() << endl;

        forAll(featurePoints, i)
        {
            label pointi = featurePoints[i];
            meshTools::writeOBJ(str, mesh.points()[pointi]);
        }
    }


    // Get all feature edges.
    labelList meshEdges
    (
        allBoundary.meshEdges
        (
            mesh.edges(),
            mesh.cellEdges(),
            SubList<label>
            (
                mesh.faceOwner(),
                allBoundary.size(),
                mesh.nInternalFaces()
            )
        )
    );

    DynamicList<label> allFeatureEdges(isFeatureEdge.size());
    forAll(isFeatureEdge, edgeI)
    {
        if (isFeatureEdge[edgeI])
        {
            // Store in mesh edge label.
            allFeatureEdges.append(meshEdges[edgeI]);
        }
    }
    featureEdges.transfer(allFeatureEdges);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyDualMesh::~polyDualMesh()
{}


// ************************************************************************* //
