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

#include "addPatchCellLayer.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyTopoChangeMap.H"
#include "syncTools.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(addPatchCellLayer, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::addPatchCellLayer::nbrFace
(
    const labelListList& edgeFaces,
    const label edgeI,
    const label facei
)
{
    const labelList& eFaces = edgeFaces[edgeI];

    if (eFaces.size() == 2)
    {
        return (eFaces[0] != facei ? eFaces[0] : eFaces[1]);
    }
    else
    {
        return -1;
    }
}


void Foam::addPatchCellLayer::addVertex
(
    const label pointi,
    face& f,
    label& fp
)
{
    if (fp == 0)
    {
        f[fp++] = pointi;
    }
    else
    {
        if (f[fp-1] != pointi && f[0] != pointi)
        {
            f[fp++] = pointi;
        }
    }
}


bool Foam::addPatchCellLayer::sameEdgeNeighbour
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const boolList& doneEdge,
    const label thisGlobalFacei,
    const label nbrGlobalFacei,
    const label edgeI
) const
{
    const edge& e = pp.edges()[edgeI];

    return
        !doneEdge[edgeI]                            // not yet handled
     && (
            addedPoints_[e[0]].size()               // is extruded
         || addedPoints_[e[1]].size()
        )
     && (
            nbrFace(globalEdgeFaces, edgeI, thisGlobalFacei)
         == nbrGlobalFacei  // is to same neighbour
        );
}


Foam::labelPair Foam::addPatchCellLayer::getEdgeString
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const boolList& doneEdge,
    const label patchFacei,
    const label globalFacei
) const
{
    const labelList& fEdges = pp.faceEdges()[patchFacei];

    label startFp = -1;
    label endFp = -1;

    // Get edge that hasn't been done yet but needs extrusion
    forAll(fEdges, fp)
    {
        label edgeI = fEdges[fp];
        const edge& e = pp.edges()[edgeI];

        if
        (
            !doneEdge[edgeI]
         && ( addedPoints_[e[0]].size() || addedPoints_[e[1]].size() )
        )
        {
            startFp = fp;
            break;
        }
    }

    if (startFp != -1)
    {
        // We found an edge that needs extruding but hasn't been done yet.
        // Now find the face on the other side
        label nbrGlobalFacei = nbrFace
        (
            globalEdgeFaces,
            fEdges[startFp],
            globalFacei
        );

        if (nbrGlobalFacei == -1)
        {
            // Proper boundary edge. Only extrude single edge.
            endFp = startFp;
        }
        else
        {
            // Search back for edge
            // - which hasn't been handled yet
            // - with same neighbour
            // - that needs extrusion
            while (true)
            {
                label prevFp = fEdges.rcIndex(startFp);

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFacei,
                        nbrGlobalFacei,
                        fEdges[prevFp]
                    )
                )
                {
                    break;
                }
                startFp = prevFp;
            }

            // Search forward for end of string
            endFp = startFp;
            while (true)
            {
                label nextFp = fEdges.fcIndex(endFp);

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFacei,
                        nbrGlobalFacei,
                        fEdges[nextFp]
                    )
                )
                {
                    break;
                }
                endFp = nextFp;
            }
        }
    }

    return labelPair(startFp, endFp);
}


Foam::label Foam::addPatchCellLayer::addSideFace
(
    const indirectPrimitivePatch& pp,
    const labelListList& addedCells,    // per pp face the new extruded cell
    const face& newFace,
    const label newPatchID,

    const label ownFacei,               // pp face that provides owner
    const label nbrFacei,
    const label layerI,                 // layer
    const label numEdgeFaces,           // number of layers for edge
    const labelList& meshFaces,         // precalculated edgeFaces
    polyTopoChange& meshMod
) const
{
    label masterFacei = -1;

    // Zone info comes from any side patch face. Otherwise -1 since we
    // don't know what to put it in - inherit from the extruded faces?
    label zoneI = -1;   // mesh_.faceZones().whichZone(meshFacei);
    bool flip = false;

    label addedFacei = -1;

    // Is patch edge external edge of indirectPrimitivePatch?
    if (nbrFacei == -1)
    {
        // External edge so external face.

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Loop over all faces connected to edge and see if we can find a face
        // that is otherPatchID

        // Get my mesh face and its zone.
        label meshFacei = pp.addressing()[ownFacei];

        forAll(meshFaces, k)
        {
            label facei = meshFaces[k];

            if
            (
                (facei != meshFacei)
             && (patches.whichPatch(facei) == newPatchID)
            )
            {
                // Found the patch face. Use it to map from
                masterFacei = facei;

                zoneI = mesh_.faceZones().whichZone(facei);
                if (zoneI != -1)
                {
                    label index = mesh_.faceZones()[zoneI].whichFace(facei);
                    flip = mesh_.faceZones()[zoneI].flipMap()[index];
                }
                break;
            }
        }

        // Determine if different number of layer on owner and neighbour side
        // (relevant only for coupled faces). See section for internal edge
        // below.

        label layerOwn;

        if (addedCells[ownFacei].size() < numEdgeFaces)
        {
            label offset = numEdgeFaces - addedCells[ownFacei].size();
            if (layerI <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layerI - offset;
            }
        }
        else
        {
            layerOwn = layerI;
        }


        // Pout<< "Added boundary face:" << newFace
        //    << " own:" << addedCells[ownFacei][layerOwn]
        //    << " patch:" << newPatchID
        //    << endl;

        addedFacei = meshMod.addFace
        (
            newFace,                    // face
            addedCells[ownFacei][layerOwn],   // owner
            -1,                         // neighbour
            masterFacei,                // master face
            false,                      // flux flip
            newPatchID,                 // patch for face
            zoneI,                      // zone for face
            flip                        // face zone flip
        );
    }
    else
    {
        // When adding side faces we need to modify neighbour and owners
        // in region where layer mesh is stopped. Determine which side
        // has max number of faces and make sure layers match closest to
        // original pp if there are different number of layers.

        label layerNbr;
        label layerOwn;

        if (addedCells[ownFacei].size() > addedCells[nbrFacei].size())
        {
            label offset =
                addedCells[ownFacei].size() - addedCells[nbrFacei].size();

            layerOwn = layerI;

            if (layerI <= offset)
            {
                layerNbr = 0;
            }
            else
            {
                layerNbr = layerI - offset;
            }
        }
        else if (addedCells[nbrFacei].size() > addedCells[ownFacei].size())
        {
            label offset =
                addedCells[nbrFacei].size() - addedCells[ownFacei].size();

            layerNbr = layerI;

            if (layerI <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layerI - offset;
            }
        }
        else
        {
            // Same number of layers on both sides.
            layerNbr = layerI;
            layerOwn = layerI;
        }

        addedFacei = meshMod.addFace
        (
            newFace,                    // face
            addedCells[ownFacei][layerOwn],   // owner
            addedCells[nbrFacei][layerNbr],   // neighbour
            -1,                         // master face
            false,                      // flux flip
            -1,                         // patch for face
            zoneI,                      // zone for face
            flip                        // face zone flip
        );

       // Pout<< "Added internal face:" << newFace
        //    << " own:" << addedCells[ownFacei][layerOwn]
        //    << " nei:" << addedCells[nbrFacei][layerNbr]
        //    << endl;
    }

    return addedFacei;
}


Foam::label Foam::addPatchCellLayer::findProcPatch
(
    const polyMesh& mesh,
    const label nbrProcID
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(mesh.globalData().processorPatches(), i)
    {
        label patchi = mesh.globalData().processorPatches()[i];

        if
        (
            refCast<const processorPolyPatch>(patches[patchi]).neighbProcNo()
         == nbrProcID
        )
        {
            return patchi;
        }
    }
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::addPatchCellLayer::addPatchCellLayer
(
    const polyMesh& mesh,
    const bool addToMesh
)
:
    mesh_(mesh),
    addToMesh_(addToMesh),
    addedPoints_(0),
    layerFaces_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::addPatchCellLayer::addedCells
(
    const polyMesh& mesh,
    const labelListList& layerFaces
)
{
    labelListList layerCells(layerFaces.size());

    forAll(layerFaces, patchFacei)
    {
        const labelList& faceLabels = layerFaces[patchFacei];

        if (faceLabels.size())
        {
            labelList& added = layerCells[patchFacei];
            added.setSize(faceLabels.size()-1);

            for (label i = 0; i < faceLabels.size()-1; i++)
            {
                added[i] = mesh.faceNeighbour()[faceLabels[i]];
            }
        }
    }
    return layerCells;
}


Foam::labelListList Foam::addPatchCellLayer::addedCells() const
{
    return addedCells(mesh_, layerFaces_);
}


Foam::labelListList Foam::addPatchCellLayer::globalEdgeFaces
(
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const indirectPrimitivePatch& pp
)
{
    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    // From mesh edge to global face labels. Non-empty sublists only for
    // pp edges.
    labelListList globalEdgeFaces(mesh.nEdges());

    const labelListList& edgeFaces = pp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];

        const labelList& eFaces = edgeFaces[edgeI];

        // Store face and processor as unique tag.
        labelList& globalEFaces = globalEdgeFaces[meshEdgeI];
        globalEFaces.setSize(eFaces.size());
        forAll(eFaces, i)
        {
            globalEFaces[i] = globalFaces.toGlobal(pp.addressing()[eFaces[i]]);
        }
    }

    // Synchronise across coupled edges.
    syncTools::syncEdgeList
    (
        mesh,
        globalEdgeFaces,
        uniqueEqOp(),
        labelList()             // null value
    );

    // Extract pp part
    return labelListList(UIndirectList<labelList>(globalEdgeFaces, meshEdges));
}


void Foam::addPatchCellLayer::calcSidePatch
(
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const labelListList& globalEdgeFaces,
    const indirectPrimitivePatch& pp,

    labelList& sidePatchID,
    label& nPatches,
    Map<label>& nbrProcToPatch,
    Map<label>& patchToNbrProc
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    sidePatchID.setSize(pp.nEdges());
    sidePatchID = -1;

    // These also get determined but not (yet) exported:
    // - whether face is created from other face or edge
    // - what zone&orientation face should have

    labelList masterFacei(pp.nEdges(), -1);

    nPatches = patches.size();

    forAll(globalEdgeFaces, edgeI)
    {
        const labelList& eGlobalFaces = globalEdgeFaces[edgeI];
        if
        (
            eGlobalFaces.size() == 2
         && pp.edgeFaces()[edgeI].size() == 1
        )
        {
            // Locally but not globally a boundary edge. Hence a coupled
            // edge. Find the patch to use if on different
            // processors.

            label f0 = eGlobalFaces[0];
            label f1 = eGlobalFaces[1];

            label otherProci = -1;
            if (globalFaces.isLocal(f0) && !globalFaces.isLocal(f1))
            {
                otherProci = globalFaces.whichProcID(f1);
            }
            else if (!globalFaces.isLocal(f0) && globalFaces.isLocal(f1))
            {
                otherProci = globalFaces.whichProcID(f0);
            }


            if (otherProci != -1)
            {
                sidePatchID[edgeI] = findProcPatch(mesh, otherProci);
                if (sidePatchID[edgeI] == -1)
                {
                    // Cannot find a patch to processor. See if already
                    // marked for addition
                    if (nbrProcToPatch.found(otherProci))
                    {
                        sidePatchID[edgeI] = nbrProcToPatch[otherProci];
                    }
                    else
                    {
                        sidePatchID[edgeI] = nPatches;
                        nbrProcToPatch.insert(otherProci, nPatches);
                        patchToNbrProc.insert(nPatches, otherProci);
                        nPatches++;
                    }
                }
            }
        }
    }



    // Determine face properties for all other boundary edges
    // ------------------------------------------------------

    const labelListList& edgeFaces = pp.edgeFaces();

    DynamicList<label> dynMeshEdgeFaces;

    forAll(edgeFaces, edgeI)
    {
        if (edgeFaces[edgeI].size() == 1 && sidePatchID[edgeI] == -1)
        {
            // Proper, uncoupled patch edge.

            label myFacei = pp.addressing()[edgeFaces[edgeI][0]];

            // Pick up any boundary face on this edge and use its properties
            label meshEdgeI = meshEdges[edgeI];
            const labelList& meshFaces = mesh.edgeFaces
            (
                meshEdgeI,
                dynMeshEdgeFaces
            );

            forAll(meshFaces, k)
            {
                label facei = meshFaces[k];

                if (facei != myFacei && !mesh.isInternalFace(facei))
                {
                    sidePatchID[edgeI] = mesh.boundaryMesh().whichPatch(facei);
                    masterFacei[edgeI] = facei;

                    break;
                }
            }
        }
    }


    // Now hopefully every boundary edge has a side patch. Check
    if (debug)
    {
        forAll(edgeFaces, edgeI)
        {
            if (edgeFaces[edgeI].size() == 1 && sidePatchID[edgeI] == -1)
            {
                const edge& e = pp.edges()[edgeI];
                WarningInFunction
                    << "Have no sidePatchID for edge " << edgeI << " points "
                    << pp.points()[pp.meshPoints()[e[0]]]
                    << pp.points()[pp.meshPoints()[e[1]]]
                    << endl;
            }
        }
    }


    // Now we have sidepatch see if we have patchface or edge to map from
    forAll(edgeFaces, edgeI)
    {
        if
        (
            edgeFaces[edgeI].size() == 1
         && sidePatchID[edgeI] != -1
         && masterFacei[edgeI] == -1
        )
        {
            // 1. Do we have a boundary face to map from

            label myFacei = pp.addressing()[edgeFaces[edgeI][0]];

            // Pick up any boundary face on this edge and use its properties
            label meshEdgeI = meshEdges[edgeI];
            const labelList& meshFaces = mesh.edgeFaces
            (
                meshEdgeI,
                dynMeshEdgeFaces
            );

            forAll(meshFaces, k)
            {
                const label facei = meshFaces[k];

                if
                (
                    facei != myFacei
                 && !mesh.isInternalFace(facei)
                 && patches.whichPatch(facei) == sidePatchID[edgeI]
                )
                {
                    sidePatchID[edgeI] = mesh.boundaryMesh().whichPatch(facei);
                    masterFacei[edgeI] = facei;
                    break;
                }
            }
        }
    }
}


void Foam::addPatchCellLayer::setRefinement
(
    const globalIndex& globalFaces,
    const labelListList& globalEdgeFaces,
    const scalarField& expansionRatio,
    const indirectPrimitivePatch& pp,
    const labelList& sidePatchID,
    const labelList& exposedPatchID,
    const labelList& nFaceLayers,
    const labelList& nPointLayers,
    const vectorField& firstLayerDisp,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "addPatchCellLayer::setRefinement : Adding up to "
            << gMax(nPointLayers)
            << " layers of cells to indirectPrimitivePatch with "
            << pp.nPoints() << " points" << endl;
    }

    if
    (
        pp.nPoints() != firstLayerDisp.size()
     || pp.nPoints() != nPointLayers.size()
     || pp.size() != nFaceLayers.size()
    )
    {
        FatalErrorInFunction
            << "Size of new points is not same as number of points used by"
            << " the face subset" << endl
            << "  patch.nPoints:" << pp.nPoints()
            << "  displacement:" << firstLayerDisp.size()
            << "  nPointLayers:" << nPointLayers.size() << nl
            << " patch.nFaces:" << pp.size()
            << "  nFaceLayers:" << nFaceLayers.size()
            << abort(FatalError);
    }

    forAll(nPointLayers, i)
    {
        if (nPointLayers[i] < 0)
        {
            FatalErrorInFunction
                << "Illegal number of layers " << nPointLayers[i]
                << " at patch point " << i << abort(FatalError);
        }
    }
    forAll(nFaceLayers, i)
    {
        if (nFaceLayers[i] < 0)
        {
            FatalErrorInFunction
                << "Illegal number of layers " << nFaceLayers[i]
                << " at patch face " << i << abort(FatalError);
        }
    }

    forAll(globalEdgeFaces, edgeI)
    {
        if (globalEdgeFaces[edgeI].size() > 2)
        {
            const edge& e = pp.edges()[edgeI];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                FatalErrorInFunction
                    << "Trying to extrude edge "
                    << e.line(pp.localPoints())
                    << " which is non-manifold (has "
                    << globalEdgeFaces[edgeI].size()
                    << " faces using it)"
                    << abort(FatalError);
            }
        }
    }


    const labelList& meshPoints = pp.meshPoints();

    // Some storage for edge-face-addressing.
    DynamicList<label> ef;

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh_.edges(), mesh_.pointEdges()));

    if (debug)
    {
        // Check synchronisation
        // ~~~~~~~~~~~~~~~~~~~~~

        {
            labelList n(mesh_.nPoints(), 0);
            UIndirectList<label>(n, meshPoints) = nPointLayers;
            syncTools::syncPointList(mesh_, n, maxEqOp<label>(), label(0));

            // Non-synced
            forAll(meshPoints, i)
            {
                label meshPointi = meshPoints[i];

                if (n[meshPointi] != nPointLayers[i])
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "On coupled point a different nLayers:"
                        << n[meshPointi] << " was specified."
                        << abort(FatalError);
                }
            }


            // Check that nPointLayers equals the max layers of connected faces
            // (or 0). Anything else makes no sense.
            labelList nFromFace(mesh_.nPoints(), 0);
            forAll(nFaceLayers, i)
            {
                const face& f = pp[i];

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    nFromFace[pointi] = max(nFromFace[pointi], nFaceLayers[i]);
                }
            }
            syncTools::syncPointList
            (
                mesh_,
                nFromFace,
                maxEqOp<label>(),
                label(0)
            );

            forAll(nPointLayers, i)
            {
                label meshPointi = meshPoints[i];

                if
                (
                    nPointLayers[i] > 0
                 && nPointLayers[i] != nFromFace[meshPointi]
                )
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "but the max nLayers of surrounding faces is:"
                        << nFromFace[meshPointi]
                        << abort(FatalError);
                }
            }
        }

        {
            pointField d(mesh_.nPoints(), vector::max);
            UIndirectList<point>(d, meshPoints) = firstLayerDisp;
            syncTools::syncPointList
            (
                mesh_,
                d,
                minEqOp<vector>(),
                vector::max
            );

            forAll(meshPoints, i)
            {
                label meshPointi = meshPoints[i];

                if (mag(d[meshPointi] - firstLayerDisp[i]) > small)
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified displacement:" << firstLayerDisp[i]
                        << endl
                        << "On coupled point a different displacement:"
                        << d[meshPointi] << " was specified."
                        << abort(FatalError);
                }
            }
        }

        // Check that edges of pp (so ones that become boundary faces)
        // connect to only one boundary face. Guarantees uniqueness of
        // patch that they go into so if this is a coupled patch both
        // sides decide the same.
        // ~~~~~~~~~~~~~~~~~~~~~~

        for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
        {
            const edge& e = pp.edges()[edgeI];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                // Edge is to become a face

                const labelList& eFaces = pp.edgeFaces()[edgeI];

                // First check: pp should be single connected.
                if (eFaces.size() != 1)
                {
                    FatalErrorInFunction
                        << "boundary-edge-to-be-extruded:"
                        << pp.points()[meshPoints[e[0]]]
                        << pp.points()[meshPoints[e[1]]]
                        << " has more than two faces using it:" << eFaces
                        << abort(FatalError);
                }

                label myFacei = pp.addressing()[eFaces[0]];

                label meshEdgeI = meshEdges[edgeI];

                // Mesh faces using edge
                const labelList& meshFaces = mesh_.edgeFaces(meshEdgeI, ef);

                // Check that there is only one patchface using edge.
                const polyBoundaryMesh& patches = mesh_.boundaryMesh();

                label bFacei = -1;

                forAll(meshFaces, i)
                {
                    label facei = meshFaces[i];

                    if (facei != myFacei)
                    {
                        if (!mesh_.isInternalFace(facei))
                        {
                            if (bFacei == -1)
                            {
                                bFacei = facei;
                            }
                            else
                            {
                                FatalErrorInFunction
                                    << "boundary-edge-to-be-extruded:"
                                    << pp.points()[meshPoints[e[0]]]
                                    << pp.points()[meshPoints[e[1]]]
                                    << " has more than two boundary faces"
                                    << " using it:"
                                    << bFacei << " fc:"
                                    << mesh_.faceCentres()[bFacei]
                                    << " patch:" << patches.whichPatch(bFacei)
                                    << " and " << facei << " fc:"
                                    << mesh_.faceCentres()[facei]
                                    << " patch:" << patches.whichPatch(facei)
                                    << abort(FatalError);
                            }
                        }
                    }
                }
            }
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Precalculated patchID for each patch face
    labelList patchID(pp.size());

    forAll(pp, patchFacei)
    {
        label meshFacei = pp.addressing()[patchFacei];

        patchID[patchFacei] = patches.whichPatch(meshFacei);
    }


    // From master point (in patch point label) to added points (in mesh point
    // label)
    addedPoints_.setSize(pp.nPoints());

    // Mark points that do not get extruded by setting size of addedPoints_ to 0
    label nTruncated = 0;

    forAll(nPointLayers, patchPointi)
    {
        if (nPointLayers[patchPointi] > 0)
        {
            addedPoints_[patchPointi].setSize(nPointLayers[patchPointi]);
        }
        else
        {
            nTruncated++;
        }
    }

    if (debug)
    {
        Pout<< "Not adding points at " << nTruncated << " out of "
            << pp.nPoints() << " points" << endl;
    }


    //
    // Create new points
    //

    // If creating new mesh: copy existing patch points
    labelList copiedPatchPoints;
    if (!addToMesh_)
    {
        copiedPatchPoints.setSize(firstLayerDisp.size());
        forAll(firstLayerDisp, patchPointi)
        {
            if (addedPoints_[patchPointi].size())
            {
                label meshPointi = meshPoints[patchPointi];
                label zoneI = mesh_.pointZones().whichZone(meshPointi);
                copiedPatchPoints[patchPointi] = meshMod.addPoint
                (
                    mesh_.points()[meshPointi],         // point
                    -1,         // master point
                    zoneI,      // zone for point
                    true        // supports a cell
                );
            }
        }
    }


    // Create points for additional layers
    forAll(firstLayerDisp, patchPointi)
    {
        if (addedPoints_[patchPointi].size())
        {
            label meshPointi = meshPoints[patchPointi];

            label zoneI = mesh_.pointZones().whichZone(meshPointi);

            point pt = mesh_.points()[meshPointi];

            vector disp = firstLayerDisp[patchPointi];

            forAll(addedPoints_[patchPointi], i)
            {
                pt += disp;

                label addedVertI = meshMod.addPoint
                (
                    pt,         // point
                    (addToMesh_ ? meshPointi : -1), // master point
                    zoneI,      // zone for point
                    true        // supports a cell
                );

                addedPoints_[patchPointi][i] = addedVertI;

                disp *= expansionRatio[patchPointi];
            }
        }
    }


    //
    // Add cells to all boundaryFaces
    //

    labelListList addedCells(pp.size());

    forAll(pp, patchFacei)
    {
        if (nFaceLayers[patchFacei] > 0)
        {
            addedCells[patchFacei].setSize(nFaceLayers[patchFacei]);

            label meshFacei = pp.addressing()[patchFacei];

            label ownZoneI = mesh_.cellZones().whichZone
            (
                mesh_.faceOwner()[meshFacei]
            );

            for (label i = 0; i < nFaceLayers[patchFacei]; i++)
            {
                // Note: add from cell (owner of patch face) or from face?
                // for now add from cell so we can map easily.
                addedCells[patchFacei][i] = meshMod.addCell
                (
                    (addToMesh_ ? mesh_.faceOwner()[meshFacei] : -1),
                    // master
                    ownZoneI        // zone for cell
                );
            }
        }
    }



    // Create faces on top of the original patch faces.
    // These faces are created from original patch faces outwards so in order
    // of increasing cell number. So orientation should be same as original
    // patch face for them to have owner<neighbour.

    layerFaces_.setSize(pp.size());

    forAll(pp.localFaces(), patchFacei)
    {
        label meshFacei = pp.addressing()[patchFacei];

        if (addedCells[patchFacei].size())
        {
            layerFaces_[patchFacei].setSize(addedCells[patchFacei].size() + 1);

            // Get duplicated vertices on the patch face.
            const face& f = pp.localFaces()[patchFacei];

            face newFace(f.size());

            forAll(addedCells[patchFacei], i)
            {
                forAll(f, fp)
                {
                    if (addedPoints_[f[fp]].empty())
                    {
                        // Keep original point
                        newFace[fp] =
                        (
                            addToMesh_
                          ? meshPoints[f[fp]]
                          : copiedPatchPoints[f[fp]]
                        );
                    }
                    else
                    {
                        // Get new outside point
                        label offset =
                            addedPoints_[f[fp]].size()
                          - addedCells[patchFacei].size();
                        newFace[fp] = addedPoints_[f[fp]][i+offset];
                    }
                }


                // Get new neighbour
                label nei;
                label patchi;
                label zoneI = -1;
                bool flip = false;


                if (i == addedCells[patchFacei].size()-1)
                {
                    // Top layer so is patch face.
                    nei = -1;
                    patchi = patchID[patchFacei];
                    zoneI = mesh_.faceZones().whichZone(meshFacei);
                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh_.faceZones()[zoneI];
                        flip = fz.flipMap()[fz.whichFace(meshFacei)];
                    }
                }
                else
                {
                    // Internal face between layer i and i+1
                    nei = addedCells[patchFacei][i+1];
                    patchi = -1;
                }


                layerFaces_[patchFacei][i+1] = meshMod.addFace
                (
                    newFace,                    // face
                    addedCells[patchFacei][i],  // owner
                    nei,                        // neighbour
                    (addToMesh_ ? meshFacei : -1), // master face
                    false,                      // flux flip
                    patchi,                     // patch for face
                    zoneI,                      // zone for face
                    flip                        // face zone flip
                );
            }
        }
    }

    //
    // Modify old patch faces to be on the inside
    //

    if (addToMesh_)
    {
        forAll(pp, patchFacei)
        {
            if (addedCells[patchFacei].size())
            {
                label meshFacei = pp.addressing()[patchFacei];

                layerFaces_[patchFacei][0] = meshFacei;

                meshMod.modifyFace
                (
                    pp[patchFacei],                 // modified face
                    meshFacei,                      // label of face
                    mesh_.faceOwner()[meshFacei],   // owner
                    addedCells[patchFacei][0],      // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    -1, // zoneI,                    // zone for face
                    false                           // face flip in zone
                );
            }
        }
    }
    else
    {
        // If creating new mesh: reverse original faces and put them
        // in the exposed patch ID.
        forAll(pp, patchFacei)
        {
            if (nFaceLayers[patchFacei] > 0)
            {
                label meshFacei = pp.addressing()[patchFacei];
                label zoneI = mesh_.faceZones().whichZone(meshFacei);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    const faceZone& fz = mesh_.faceZones()[zoneI];
                    zoneFlip = !fz.flipMap()[fz.whichFace(meshFacei)];
                }

                // Reverse and renumber old patch face.
                face f(pp.localFaces()[patchFacei].reverseFace());
                forAll(f, fp)
                {
                    f[fp] = copiedPatchPoints[f[fp]];
                }

                layerFaces_[patchFacei][0] = meshMod.addFace
                (
                    f,                          // modified face
                    addedCells[patchFacei][0],  // owner
                    -1,                         // neighbour
                    -1,                         // masterFace
                    true,                       // face flip
                    exposedPatchID[patchFacei], // patch for face
                    zoneI,                      // zone for face
                    zoneFlip                    // face flip in zone
                );
            }
        }
    }



    //
    // Create 'side' faces, one per edge that is being extended.
    //

    const labelListList& faceEdges = pp.faceEdges();
    const faceList& localFaces = pp.localFaces();
    const edgeList& edges = pp.edges();

    // Get number of layers per edge. This is 0 if edge is not extruded;
    // max of connected faces otherwise.
    labelList edgeLayers(pp.nEdges());

    {
        // Use list over mesh.nEdges() since syncTools does not yet support
        // partial list synchronisation.
        labelList meshEdgeLayers(mesh_.nEdges(), -1);

        forAll(meshEdges, edgeI)
        {
            const edge& e = edges[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            if ((nPointLayers[e[0]] == 0) && (nPointLayers[e[1]] == 0))
            {
                meshEdgeLayers[meshEdgeI] = 0;
            }
            else
            {
                const labelList& eFaces = pp.edgeFaces()[edgeI];

                forAll(eFaces, i)
                {
                    meshEdgeLayers[meshEdgeI] = max
                    (
                        nFaceLayers[eFaces[i]],
                        meshEdgeLayers[meshEdgeI]
                    );
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            meshEdgeLayers,
            maxEqOp<label>(),
            label(0)            // initial value
        );

        forAll(meshEdges, edgeI)
        {
            edgeLayers[edgeI] = meshEdgeLayers[meshEdges[edgeI]];
        }
    }


    // Mark off which edges have been extruded
    boolList doneEdge(pp.nEdges(), false);


    // Create faces. Per face walk connected edges and find string of edges
    // between the same two faces and extrude string into a single face.
    forAll(pp, patchFacei)
    {
        const labelList& fEdges = faceEdges[patchFacei];

        forAll(fEdges, fp)
        {
            // Get string of edges that needs to be extruded as a single face.
            // Returned as indices in fEdges.
            labelPair indexPair
            (
                getEdgeString
                (
                    pp,
                    globalEdgeFaces,
                    doneEdge,
                    patchFacei,
                    globalFaces.toGlobal(pp.addressing()[patchFacei])
                )
            );

            // Pout<< "Found unextruded edges in edges:" << fEdges
            //    << " start:" << indexPair[0]
            //    << " end:" << indexPair[1]
            //    << endl;

            const label startFp = indexPair[0];
            const label endFp = indexPair[1];

            if (startFp != -1)
            {
                // Extrude edges from indexPair[0] up to indexPair[1]
                // (note indexPair = indices of edges. There is one more vertex
                //  than edges)
                const face& f = localFaces[patchFacei];

                labelList stringedVerts;
                if (endFp >= startFp)
                {
                    stringedVerts.setSize(endFp-startFp+2);
                }
                else
                {
                    stringedVerts.setSize(endFp+f.size()-startFp+2);
                }

                label fp = startFp;

                for (label i = 0; i < stringedVerts.size()-1; i++)
                {
                    stringedVerts[i] = f[fp];
                    doneEdge[fEdges[fp]] = true;
                    fp = f.fcIndex(fp);
                }
                stringedVerts.last() = f[fp];


                // Now stringedVerts contains the vertices in order of face f.
                // This is consistent with the order if f becomes the owner cell
                // and nbrFacei the neighbour cell. Note that the cells get
                // added in order of pp so we can just use face ordering and
                // because we loop in incrementing order as well we will
                // always have nbrFacei > patchFacei.

                label startEdgeI = fEdges[startFp];

                label meshEdgeI = meshEdges[startEdgeI];

                label numEdgeSideFaces = edgeLayers[startEdgeI];

                for (label i = 0; i < numEdgeSideFaces; i++)
                {
                    label vEnd = stringedVerts.last();
                    label vStart = stringedVerts[0];

                    // calculate number of points making up a face
                    label newFp = 2*stringedVerts.size();

                    if (i == 0)
                    {
                        // layer 0 gets all the truncation of neighbouring
                        // faces with more layers.
                        if (addedPoints_[vEnd].size())
                        {
                            newFp +=
                                addedPoints_[vEnd].size() - numEdgeSideFaces;
                        }
                        if (addedPoints_[vStart].size())
                        {
                            newFp +=
                                addedPoints_[vStart].size() - numEdgeSideFaces;
                        }
                    }

                    face newFace(newFp);

                    newFp = 0;

                    // For layer 0 get pp points, for all other layers get
                    // points of layer-1.
                    if (i == 0)
                    {
                        forAll(stringedVerts, stringedI)
                        {
                            label v = stringedVerts[stringedI];
                            addVertex
                            (
                                (
                                    addToMesh_
                                  ? meshPoints[v]
                                  : copiedPatchPoints[v]
                                ),
                                newFace,
                                newFp
                            );
                        }
                    }
                    else
                    {
                        forAll(stringedVerts, stringedI)
                        {
                            label v = stringedVerts[stringedI];
                            if (addedPoints_[v].size())
                            {
                                label offset =
                                    addedPoints_[v].size() - numEdgeSideFaces;
                                addVertex
                                (
                                    addedPoints_[v][i+offset-1],
                                    newFace,
                                    newFp
                                );
                            }
                            else
                            {
                                addVertex
                                (
                                    (
                                        addToMesh_
                                      ? meshPoints[v]
                                      : copiedPatchPoints[v]
                                    ),
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    // add points between stringed vertices (end)
                    if (numEdgeSideFaces < addedPoints_[vEnd].size())
                    {
                        if (i == 0 && addedPoints_[vEnd].size())
                        {
                            label offset =
                                addedPoints_[vEnd].size() - numEdgeSideFaces;
                            for (label ioff = 0; ioff < offset; ioff++)
                            {
                                addVertex
                                (
                                    addedPoints_[vEnd][ioff],
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    forAllReverse(stringedVerts, stringedI)
                    {
                        label v = stringedVerts[stringedI];
                        if (addedPoints_[v].size())
                        {
                            label offset =
                                addedPoints_[v].size() - numEdgeSideFaces;
                            addVertex
                            (
                                addedPoints_[v][i+offset],
                                newFace,
                                newFp
                            );
                        }
                        else
                        {
                            addVertex
                            (
                                (
                                    addToMesh_
                                  ? meshPoints[v]
                                  : copiedPatchPoints[v]
                                ),
                                newFace,
                                newFp
                            );
                        }
                    }


                    // add points between stringed vertices (start)
                    if (numEdgeSideFaces < addedPoints_[vStart].size())
                    {
                        if (i == 0 && addedPoints_[vStart].size())
                        {
                            label offset =
                                addedPoints_[vStart].size() - numEdgeSideFaces;
                            for (label ioff = offset-1; ioff >= 0; ioff--)
                            {
                                addVertex
                                (
                                    addedPoints_[vStart][ioff],
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    if (newFp >= 3)
                    {
                        // Add face in between faces patchFacei and nbrFacei
                        // (possibly -1 for external edges)

                        newFace.setSize(newFp);

                        if (debug)
                        {
                            labelHashSet verts(2*newFace.size());
                            forAll(newFace, fp)
                            {
                                if (!verts.insert(newFace[fp]))
                                {
                                    FatalErrorInFunction
                                        << "Duplicate vertex in face"
                                        << " to be added." << nl
                                        << "newFace:" << newFace << nl
                                        << "points:"
                                        <<  UIndirectList<point>
                                            (
                                                meshMod.points(),
                                                newFace
                                            ) << nl
                                        << "Layer:" << i
                                        << " out of:" << numEdgeSideFaces << nl
                                        << "ExtrudeEdge:" << meshEdgeI
                                        << " at:"
                                        <<  mesh_.edges()[meshEdgeI].line
                                            (
                                                mesh_.points()
                                            ) << nl
                                        << "string:" << stringedVerts
                                        << "stringpoints:"
                                        << UIndirectList<point>
                                            (
                                                pp.localPoints(),
                                                stringedVerts
                                            ) << nl
                                        << "stringNLayers:"
                                        <<  UIndirectList<label>
                                            (
                                                nPointLayers,
                                                stringedVerts
                                            ) << nl
                                        << abort(FatalError);
                                }
                            }
                        }

                        label nbrFacei = nbrFace
                        (
                            pp.edgeFaces(),
                            startEdgeI,
                            patchFacei
                        );

                        const labelList& meshFaces = mesh_.edgeFaces
                        (
                            meshEdgeI,
                            ef
                        );

                        addSideFace
                        (
                            pp,
                            addedCells,

                            newFace,                // vertices of new face
                            sidePatchID[startEdgeI],// -1 or patch for face

                            patchFacei,
                            nbrFacei,
                            i,                  // layer
                            numEdgeSideFaces,   // num layers
                            meshFaces,          // edgeFaces
                            meshMod
                        );
                    }
                }
            }
        }
    }
}


void Foam::addPatchCellLayer::topoChange
(
    const polyTopoChangeMap& map,
    const labelList& faceMap,   // new to old patch faces
    const labelList& pointMap   // new to old patch points
)
{
    {
        labelListList newAddedPoints(pointMap.size());

        forAll(newAddedPoints, newPointi)
        {
            label oldPointi = pointMap[newPointi];

            const labelList& added = addedPoints_[oldPointi];

            labelList& newAdded = newAddedPoints[newPointi];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newPointi = map.reversePointMap()[added[i]];

                if (newPointi >= 0)
                {
                    newAdded[newI++] = newPointi;
                }
            }
            newAdded.setSize(newI);
        }
        addedPoints_.transfer(newAddedPoints);
    }

    {
        labelListList newLayerFaces(faceMap.size());

        forAll(newLayerFaces, newFacei)
        {
            label oldFacei = faceMap[newFacei];

            const labelList& added = layerFaces_[oldFacei];

            labelList& newAdded = newLayerFaces[newFacei];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newFacei = map.reverseFaceMap()[added[i]];

                if (newFacei >= 0)
                {
                    newAdded[newI++] = newFacei;
                }
            }
            newAdded.setSize(newI);
        }
        layerFaces_.transfer(newLayerFaces);
    }
}


// ************************************************************************* //
