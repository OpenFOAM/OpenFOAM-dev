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

\*---------------------------------------------------------------------------*/

#include "BiIndirectList.H"
#include "removePoints.H"
#include "PstreamReduceOps.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyRemovePoint.H"
#include "polyAddPoint.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "faceSet.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(removePoints, 0);

// Combine-reduce operator to combine data on faces. Takes care
// of reverse orientation on coupled face.
template<class T, template<class> class CombineOp>
class faceEqOp
{

public:

    void operator()(List<T>& x, const List<T>& y) const
    {
        if (y.size() > 0)
        {
            if (x.size() == 0)
            {
                x = y;
            }
            else
            {
                label j = 0;
                forAll(x, i)
                {
                    CombineOp<T>()(x[i], y[j]);
                    j = y.rcIndex(j);
                }
            }
        }
    }
};

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::removePoints::modifyFace
(
    const label facei,
    const face& newFace,
    polyTopoChange& meshMod
) const
{
    // Get other face data.
    label patchi = -1;
    label owner = mesh_.faceOwner()[facei];
    label neighbour = -1;

    if (mesh_.isInternalFace(facei))
    {
        neighbour = mesh_.faceNeighbour()[facei];
    }
    else
    {
        patchi = mesh_.boundaryMesh().whichPatch(facei);
    }

    label zoneID = mesh_.faceZones().whichZone(facei);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            newFace,        // modified face
            facei,          // label of face being modified
            owner,          // owner
            neighbour,      // neighbour
            false,          // face flip
            patchi,         // patch for face
            false,          // remove from zone
            zoneID,         // zone for face
            zoneFlip        // face flip in zone
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::removePoints::removePoints
(
    const polyMesh& mesh,
    const bool undoable
)
:
    mesh_(mesh),
    undoable_(undoable)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::removePoints::countPointUsage
(
    const scalar minCos,
    boolList& pointCanBeDeleted
) const
{
    // Containers to store two edges per point:
    // -1   : not filled
    // >= 0 : edge label
    // -2   : more than two edges for point
    labelList edge0(mesh_.nPoints(), -1);
    labelList edge1(mesh_.nPoints(), -1);

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        forAll(e, eI)
        {
            label pointi = e[eI];

            if (edge0[pointi] == -2)
            {
                // Already too many edges
            }
            else if (edge0[pointi] == -1)
            {
                // Store first edge using point
                edge0[pointi] = edgeI;
            }
            else
            {
                // Already one edge using point. Check second container.

                if (edge1[pointi] == -1)
                {
                    // Store second edge using point
                    edge1[pointi] = edgeI;
                }
                else
                {
                    // Third edge using point. Mark.
                    edge0[pointi] = -2;
                    edge1[pointi] = -2;
                }
            }
        }
    }


    // Check the ones used by only 2 edges that these are sufficiently in line.
    const pointField& points = mesh_.points();

    pointCanBeDeleted.setSize(mesh_.nPoints());
    pointCanBeDeleted = false;
    //label nDeleted = 0;

    forAll(edge0, pointi)
    {
        if (edge0[pointi] >= 0 && edge1[pointi] >= 0)
        {
            // Point used by two edges exactly

            const edge& e0 = edges[edge0[pointi]];
            const edge& e1 = edges[edge1[pointi]];

            label common = e0.commonVertex(e1);
            label vLeft = e0.otherVertex(common);
            label vRight = e1.otherVertex(common);

            vector e0Vec = points[common] - points[vLeft];
            e0Vec /= mag(e0Vec) + vSmall;

            vector e1Vec = points[vRight] - points[common];
            e1Vec /= mag(e1Vec) + vSmall;

            if ((e0Vec & e1Vec) > minCos)
            {
                pointCanBeDeleted[pointi] = true;
                //nDeleted++;
            }
        }
        else if (edge0[pointi] == -1)
        {
            // point not used at all
            pointCanBeDeleted[pointi] = true;
            //nDeleted++;
        }
    }
    edge0.clear();
    edge1.clear();


    // Protect any points on faces that would collapse down to nothing
    // No particular intelligence so might protect too many points
    forAll(mesh_.faces(), facei)
    {
        const face& f = mesh_.faces()[facei];

        label nCollapse = 0;
        forAll(f, fp)
        {
            if (pointCanBeDeleted[f[fp]])
            {
                nCollapse++;
            }
        }

        if ((f.size() - nCollapse) < 3)
        {
            // Just unmark enough points
            forAll(f, fp)
            {
                if (pointCanBeDeleted[f[fp]])
                {
                    pointCanBeDeleted[f[fp]] = false;
                    --nCollapse;
                    if (nCollapse == 0)
                    {
                        break;
                    }
                }
            }
        }
    }


    // Point can be deleted only if all processors want to delete it
    syncTools::syncPointList
    (
        mesh_,
        pointCanBeDeleted,
        andEqOp<bool>(),
        true                // null value
    );

    label nDeleted = 0;
    forAll(pointCanBeDeleted, pointi)
    {
        if (pointCanBeDeleted[pointi])
        {
            nDeleted++;
        }
    }

    return returnReduce(nDeleted, sumOp<label>());
}


void Foam::removePoints::setRefinement
(
    const boolList& pointCanBeDeleted,
    polyTopoChange& meshMod
)
{
    // Count deleted points
    label nDeleted = 0;
    forAll(pointCanBeDeleted, pointi)
    {
        if (pointCanBeDeleted[pointi])
        {
            nDeleted++;
        }
    }

    // Faces (in mesh face labels) affected by points removed. Will hopefully
    // be only a few.
    labelHashSet facesAffected(4*nDeleted);


    // Undo: from global mesh point to index in savedPoint_
    Map<label> pointToSaved;

    // Size undo storage
    if (undoable_)
    {
        savedPoints_.setSize(nDeleted);
        pointToSaved.resize(2*nDeleted);
    }


    // Remove points
    // ~~~~~~~~~~~~~

    nDeleted = 0;

    forAll(pointCanBeDeleted, pointi)
    {
        if (pointCanBeDeleted[pointi])
        {
            if (undoable_)
            {
                pointToSaved.insert(pointi, nDeleted);
                savedPoints_[nDeleted++] = mesh_.points()[pointi];
            }
            meshMod.setAction(polyRemovePoint(pointi));

            // Store faces affected
            const labelList& pFaces = mesh_.pointFaces()[pointi];

            forAll(pFaces, i)
            {
                facesAffected.insert(pFaces[i]);
            }
        }
    }



    // Update faces
    // ~~~~~~~~~~~~


    if (undoable_)
    {
        savedFaceLabels_.setSize(facesAffected.size());
        savedFaces_.setSize(facesAffected.size());
    }
    label nSaved = 0;

    forAllConstIter(labelHashSet, facesAffected, iter)
    {
        label facei = iter.key();

        const face& f = mesh_.faces()[facei];

        face newFace(f.size());

        label newI = 0;

        forAll(f, fp)
        {
            label pointi = f[fp];

            if (!pointCanBeDeleted[pointi])
            {
                newFace[newI++] = pointi;
            }
        }
        newFace.setSize(newI);

        // Actually change the face to the new vertices
        modifyFace(facei, newFace, meshMod);

        // Save the face. Negative indices are into savedPoints_
        if (undoable_)
        {
            savedFaceLabels_[nSaved] = facei;

            face& savedFace = savedFaces_[nSaved++];
            savedFace.setSize(f.size());

            forAll(f, fp)
            {
                label pointi = f[fp];

                if (pointCanBeDeleted[pointi])
                {
                    savedFace[fp] = -pointToSaved[pointi]-1;
                }
                else
                {
                    savedFace[fp] = pointi;
                }
            }
        }
    }

    if (undoable_)
    {
        // DEBUG: Compare the stored faces with the current ones.
        if (debug)
        {
            forAll(savedFaceLabels_, saveI)
            {
                // Points from the mesh
                List<point> meshPoints
                (
                    UIndirectList<point>
                    (
                        mesh_.points(),
                        mesh_.faces()[savedFaceLabels_[saveI]]  // mesh face
                    )
                );

                // Points from the saved face
                List<point> keptPoints
                (
                    BiIndirectList<point>
                    (
                        mesh_.points(),
                        savedPoints_,
                        savedFaces_[saveI]  // saved face
                    )
                );

                if (meshPoints != keptPoints)
                {
                    FatalErrorInFunction
                        << "facei:" << savedFaceLabels_[saveI] << nl
                        << "meshPoints:" << meshPoints << nl
                        << "keptPoints:" << keptPoints << nl
                        << abort(FatalError);
                }
            }
        }
    }
}


void Foam::removePoints::updateMesh(const mapPolyMesh& map)
{
    if (undoable_)
    {
        forAll(savedFaceLabels_, localI)
        {
            if (savedFaceLabels_[localI] >= 0)
            {
                label newFacei = map.reverseFaceMap()[savedFaceLabels_[localI]];

                if (newFacei == -1)
                {
                    FatalErrorInFunction
                        << "Old face " << savedFaceLabels_[localI]
                        << " seems to have disappeared."
                        << abort(FatalError);
                }
                savedFaceLabels_[localI] = newFacei;
            }
        }


        // Renumber mesh vertices (indices >=0). Leave saved vertices
        // (<0) intact.
        forAll(savedFaces_, i)
        {
            face& f = savedFaces_[i];

            forAll(f, fp)
            {
                label pointi = f[fp];

                if (pointi >= 0)
                {
                    f[fp] = map.reversePointMap()[pointi];

                    if (f[fp] == -1)
                    {
                        FatalErrorInFunction
                            << "Old point " << pointi
                            << " seems to have disappeared."
                            << abort(FatalError);
                    }
                }
            }
        }


        // DEBUG: Compare the stored faces with the current ones.
        if (debug)
        {
            forAll(savedFaceLabels_, saveI)
            {
                if (savedFaceLabels_[saveI] >= 0)
                {
                    const face& f = mesh_.faces()[savedFaceLabels_[saveI]];

                    // Get kept points of saved faces.
                    const face& savedFace = savedFaces_[saveI];

                    face keptFace(savedFace.size());
                    label keptFp = 0;

                    forAll(savedFace, fp)
                    {
                        label pointi = savedFace[fp];

                        if (pointi >= 0)
                        {
                            keptFace[keptFp++] = pointi;
                        }
                    }
                    keptFace.setSize(keptFp);

                    // Compare as faces (since f might have rotated and
                    // face::operator== takes care of this)
                    if (keptFace != f)
                    {
                        FatalErrorInFunction
                            << "facei:" << savedFaceLabels_[saveI] << nl
                            << "face:" << f << nl
                            << "keptFace:" << keptFace << nl
                            << "saved points:"
                            <<  BiIndirectList<point>
                                (
                                    mesh_.points(),
                                    savedPoints_,
                                    savedFace
                                )() << nl
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


// Given list of faces to undo picks up the local indices of the faces
// to restore. Additionally it also picks up all the faces that use
// any of the deleted points.
void Foam::removePoints::getUnrefimentSet
(
    const labelList& undoFaces,
    labelList& localFaces,
    labelList& localPoints
) const
{
    if (!undoable_)
    {
        FatalErrorInFunction
            << "removePoints not constructed with"
            << " unrefinement capability."
            << abort(FatalError);
    }

    if (debug)
    {
        // Check if synced.
        faceSet undoFacesSet(mesh_, "undoFacesSet", undoFaces);
        label sz = undoFacesSet.size();

        undoFacesSet.sync(mesh_);
        if (sz != undoFacesSet.size())
        {
            FatalErrorInFunction
                << "undoFaces not synchronised across coupled faces." << endl
                << "Size before sync:" << sz
                << "  after sync:" << undoFacesSet.size()
                << abort(FatalError);
        }
    }


    // Problem: input undoFaces are synced. But e.g.
    // two faces, A (uncoupled) and B(coupled) share a deleted point. A gets
    // passed in to be restored. Now picking up the deleted point and the
    // original faces using it picks up face B. But not its coupled neighbour!
    // Problem is that we cannot easily synchronise the deleted points
    // (e.g. syncPointList) since they're not in the mesh anymore - only the
    // faces are. So we have to collect the points-to-restore as indices
    // in the faces (which is information we can synchronise)



    // Mark points-to-restore
    labelHashSet localPointsSet(undoFaces.size());

    {
        // Create set of faces to be restored
        labelHashSet undoFacesSet(undoFaces);

        forAll(savedFaceLabels_, saveI)
        {
            if (savedFaceLabels_[saveI] < 0)
            {
                FatalErrorInFunction
                    << "Illegal face label " << savedFaceLabels_[saveI]
                    << " at index " << saveI
                    << abort(FatalError);
            }

            if (undoFacesSet.found(savedFaceLabels_[saveI]))
            {
                const face& savedFace = savedFaces_[saveI];

                forAll(savedFace, fp)
                {
                    if (savedFace[fp] < 0)
                    {
                        label savedPointi = -savedFace[fp]-1;

                        if (savedPoints_[savedPointi] == vector::max)
                        {
                            FatalErrorInFunction
                                << "Trying to restore point " << savedPointi
                                << " from mesh face " << savedFaceLabels_[saveI]
                                << " saved face:" << savedFace
                                << " which has already been undone."
                                << abort(FatalError);
                        }

                        localPointsSet.insert(savedPointi);
                    }
                }
            }
        }


        // Per boundary face, per index in face whether the point needs
        // restoring. Note that this is over all saved faces, not just over
        // the ones in undoFaces.

        boolListList faceVertexRestore(mesh_.nFaces()-mesh_.nInternalFaces());

        // Populate with my local points-to-restore.
        forAll(savedFaces_, saveI)
        {
            label bFacei = savedFaceLabels_[saveI] - mesh_.nInternalFaces();

            if (bFacei >= 0)
            {
                const face& savedFace = savedFaces_[saveI];

                boolList& fRestore = faceVertexRestore[bFacei];

                fRestore.setSize(savedFace.size());
                fRestore = false;

                forAll(savedFace, fp)
                {
                    if (savedFace[fp] < 0)
                    {
                        label savedPointi = -savedFace[fp]-1;

                        if (localPointsSet.found(savedPointi))
                        {
                            fRestore[fp] = true;
                        }
                    }
                }
            }
        }

        syncTools::syncBoundaryFaceList
        (
            mesh_,
            faceVertexRestore,
            faceEqOp<bool, orEqOp>(),   // special operator to handle faces
            Foam::dummyTransform()      // no transformation
        );

        // So now if any of the points-to-restore is used by any coupled face
        // anywhere the corresponding index in faceVertexRestore will be set.

        // Now combine the localPointSet and the (sychronised)
        // boundary-points-to-restore.

        forAll(savedFaces_, saveI)
        {
            label bFacei = savedFaceLabels_[saveI] - mesh_.nInternalFaces();

            if (bFacei >= 0)
            {
                const boolList& fRestore = faceVertexRestore[bFacei];

                const face& savedFace = savedFaces_[saveI];

                forAll(fRestore, fp)
                {
                    // Does neighbour require point restored?
                    if (fRestore[fp])
                    {
                        if (savedFace[fp] >= 0)
                        {
                            FatalErrorInFunction
                                << "Problem: on coupled face:"
                                << savedFaceLabels_[saveI]
                                << " fc:"
                                << mesh_.faceCentres()[savedFaceLabels_[saveI]]
                                << endl
                                << " my neighbour tries to restore the vertex"
                                << " at index " << fp
                                << " whereas my saved face:" << savedFace
                                << " does not indicate a deleted vertex"
                                << " at that position."
                                << abort(FatalError);
                        }

                        label savedPointi = -savedFace[fp]-1;

                        localPointsSet.insert(savedPointi);
                    }
                }
            }
        }
    }

    localPoints = localPointsSet.toc();


    // Collect all saved faces using any localPointsSet
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet localFacesSet(2*undoFaces.size());

    forAll(savedFaces_, saveI)
    {
        const face& savedFace = savedFaces_[saveI];

        forAll(savedFace, fp)
        {
            if (savedFace[fp] < 0)
            {
                label savedPointi = -savedFace[fp]-1;

                if (localPointsSet.found(savedPointi))
                {
                    localFacesSet.insert(saveI);
                }
            }
        }
    }
    localFaces = localFacesSet.toc();


    // Note that at this point the localFaces to restore can contain points
    // that are not going to be restored! The localFaces though will
    // be guaranteed to be all the ones affected by the restoration of the
    // localPoints.
}


void Foam::removePoints::setUnrefinement
(
    const labelList& localFaces,
    const labelList& localPoints,
    polyTopoChange& meshMod
)
{
    if (!undoable_)
    {
        FatalErrorInFunction
            << "removePoints not constructed with"
            << " unrefinement capability."
            << abort(FatalError);
    }


    // Per savedPoint -1 or the restored point label
    labelList addedPoints(savedPoints_.size(), -1);

    forAll(localPoints, i)
    {
        label localI = localPoints[i];

        if (savedPoints_[localI] == vector::max)
        {
            FatalErrorInFunction
                << "Saved point " << localI << " already restored!"
                << abort(FatalError);
        }

        addedPoints[localI] = meshMod.setAction
        (
            polyAddPoint
            (
                savedPoints_[localI],   // point
                -1,                     // master point
                -1,                     // zone for point
                true                    // supports a cell
            )
        );

        // Mark the restored points so they are not restored again.
        savedPoints_[localI] = vector::max;
    }

    forAll(localFaces, i)
    {
        label saveI = localFaces[i];

        // Modify indices into saved points (so < 0) to point to the
        // added points.
        face& savedFace = savedFaces_[saveI];

        face newFace(savedFace.size());
        label newFp = 0;

        bool hasSavedPoints = false;

        forAll(savedFace, fp)
        {
            if (savedFace[fp] < 0)
            {
                label addedPointi = addedPoints[-savedFace[fp]-1];

                if (addedPointi != -1)
                {
                    savedFace[fp] = addedPointi;
                    newFace[newFp++] = addedPointi;
                }
                else
                {
                    hasSavedPoints = true;
                }
            }
            else
            {
                newFace[newFp++] = savedFace[fp];
            }
        }
        newFace.setSize(newFp);

        modifyFace(savedFaceLabels_[saveI], newFace, meshMod);

        if (!hasSavedPoints)
        {
            // Face fully restored. Mark for compaction later on
            savedFaceLabels_[saveI] = -1;
            savedFaces_[saveI].clear();
        }
    }


    // Compact face labels
    label newSaveI = 0;

    forAll(savedFaceLabels_, saveI)
    {
        if (savedFaceLabels_[saveI] != -1)
        {
            if (newSaveI != saveI)
            {
                savedFaceLabels_[newSaveI] = savedFaceLabels_[saveI];
                savedFaces_[newSaveI].transfer(savedFaces_[saveI]);
            }
            newSaveI++;
        }
    }

    savedFaceLabels_.setSize(newSaveI);
    savedFaces_.setSize(newSaveI);


    // Check that all faces have been restored that use any restored points
    if (debug)
    {
        forAll(savedFaceLabels_, saveI)
        {
            const face& savedFace = savedFaces_[saveI];

            forAll(savedFace, fp)
            {
                if (savedFace[fp] < 0)
                {
                    label addedPointi = addedPoints[-savedFace[fp]-1];

                    if (addedPointi != -1)
                    {
                        FatalErrorInFunction
                            << "Face:" << savedFaceLabels_[saveI]
                            << " savedVerts:" << savedFace
                            << " uses restored point:" << -savedFace[fp]-1
                            << " with new pointlabel:" << addedPointi
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
