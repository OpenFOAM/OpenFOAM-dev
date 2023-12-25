/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "mergePatchPairs.H"
#include "fvMesh.H"
#include "polyPatchIntersection.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitch(mergePatchPairs, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::mergePatchPairs::findPatchIndex(const word& patchName) const
{
    label patchIndex = mesh_.boundaryMesh().findIndex(patchName);

    if (patchIndex == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << patchName << exit(FatalError);
    }

    return patchIndex;
}


Foam::Pair<Foam::label> Foam::mergePatchPairs::findPatchIndices
(
    const Pair<word>& patchNamePair
) const
{
    return Pair<label>
    {
        findPatchIndex(patchNamePair.first()),
        findPatchIndex(patchNamePair.second())
    };
}


Foam::List<Foam::Pair<Foam::label>> Foam::mergePatchPairs::patchPairs
(
    const List<Pair<word>>& patchNamePairs
) const
{
    List<Pair<label>> patchPairs(patchNamePairs.size());

    forAll(patchNamePairs, pairi)
    {
        patchPairs[pairi] = findPatchIndices(patchNamePairs[pairi]);
    }

    return patchPairs;
}


void Foam::mergePatchPairs::removePoints
(
    polyTopoChange& meshMod,
    const List<Pair<label>>& patchPairs
) const
{
    boolList removedPoints(mesh_.nPoints(), false);

    forAll(patchPairs, ppi)
    {
        UIndirectList<bool>
        (
            removedPoints,
            mesh_.boundaryMesh()[patchPairs[ppi].first()].meshPoints()
        ) = true;

        UIndirectList<bool>
        (
            removedPoints,
            mesh_.boundaryMesh()[patchPairs[ppi].second()].meshPoints()
        ) = true;
    }

    forAll(removedPoints, pi)
    {
        if (removedPoints[pi])
        {
            meshMod.removePoint(pi, -1);
        }
    }
}


void Foam::mergePatchPairs::addPoints
(
    polyTopoChange& meshMod,
    const polyPatchIntersection& intersection
) const
{
    const pointField& intersectionPoints = intersection.points();
    newPoints_.setSize(intersectionPoints.size(), -1);

    forAll(intersectionPoints, pi)
    {
        const label pointi = meshMod.addPoint
        (
            intersectionPoints[pi],
            -1,
            -1,
            false
        );

        if (debug)
        {
            Info<< "Adding point "
                << pointi << " "
                << intersectionPoints[pi]
                << endl;
        }

        newPoints_[pi] = pointi;
    }
}


void Foam::mergePatchPairs::removeFaces
(
    polyTopoChange& meshMod,
    const polyPatchIntersection& intersection
) const
{
    const polyPatch& srcPatch = intersection.srcPatch();
    const polyPatch& tgtPatch = intersection.tgtPatch();

    // Remove existing source face
    forAll(srcPatch, fi)
    {
        meshMod.removeFace(srcPatch.start() + fi, -1);
    }

    // Remove existing target face
    forAll(tgtPatch, fi)
    {
        meshMod.removeFace(tgtPatch.start() + fi, -1);
    }
}


Foam::face Foam::mergePatchPairs::mapFace(const face& f) const
{
    face newFace(f);

    forAll(newFace, fpi)
    {
        newFace[fpi] = newPoints_[newFace[fpi]];
    }

    return newFace;
}


void Foam::mergePatchPairs::addFaces
(
    polyTopoChange& meshMod,
    const polyPatchIntersection& intersection
) const
{
    const faceList& faces = intersection.faces();
    const DynamicList<label>& faceSrcFaces = intersection.faceSrcFaces();
    const DynamicList<label>& faceTgtFaces = intersection.faceTgtFaces();

    const label srcPatchi = intersection.srcPatch().index();
    const label srcPatchStart = intersection.srcPatch().start();

    const label tgtPatchi = intersection.tgtPatch().index();
    const label tgtPatchStart = intersection.tgtPatch().start();

    forAll(faces, fi)
    {
        const face f = mapFace(faces[fi]);

        if (faceSrcFaces[fi] != -1 && faceTgtFaces[fi] != -1)
        {
            const label srcFacei = srcPatchStart + faceSrcFaces[fi];
            const label tgtFacei = tgtPatchStart + faceTgtFaces[fi];

            const label srcOwn = mesh_.faceOwner()[srcFacei];
            const label tgtOwn = mesh_.faceOwner()[tgtFacei];

            if (srcOwn < tgtOwn)
            {
                meshMod.addFace
                (
                    f,          // Face to add
                    srcOwn,     // Owner cell
                    tgtOwn,     // Neighbour cell
                    -1,         // Master point index
                    -1,         // Master edge index
                    srcFacei,   // Master face index
                    false,      // Flip
                    -1,         // Patch index
                    -1,         // Zone index
                    false       // Zone sign
                );

                if (debug)
                {
                    Info<< "Adding internal face " << f
                        << " owner celli:" << srcOwn
                        << " neighbour celli:" << tgtOwn << endl;
                }
            }
            else
            {
                meshMod.addFace
                (
                    f.reverseFace(), // Face to add
                    tgtOwn,     // Owner cell
                    srcOwn,     // Neighbour cell
                    -1,         // Master point index
                    -1,         // Master edge index
                    tgtFacei,   // Master face index
                    false,      // Flip
                    -1,         // Patch index
                    -1,         // Zone index
                    false       // Zone sign
                );

                if (debug)
                {
                    Info<< "Adding internal face " << f
                        << " owner celli:" << tgtOwn
                        << " neighbour celli:" << srcOwn << endl;
                }
            }
        }
        else if (faceSrcFaces[fi] != -1 && faceTgtFaces[fi] == -1)
        {
            const label srcFacei = srcPatchStart + faceSrcFaces[fi];
            const label srcOwn = mesh_.faceOwner()[srcFacei];

            // Get source face zone info
            const label zoneIndex = mesh_.faceZones().whichZone(srcFacei);
            bool zoneFlip = false;
            if (zoneIndex >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneIndex];
                zoneFlip = fZone.flipMap()[fZone.whichFace(srcFacei)];
            }

            meshMod.addFace
            (
                f,          // Face to add
                srcOwn,     // Owner cell
                -1,         // Neighbour cell
                -1,         // Master point index
                -1,         // Master edge index
                srcFacei,   // Master face index
                false,      // Flip
                srcPatchi,  // Patch index
                zoneIndex,     // Zone index
                zoneFlip    // Zone sign
            );

            if (debug)
            {
                Info<< "Adding patch face " << f
                    << " owner celli:" << srcOwn
                    << " patchi:" << srcPatchi << endl;
            }
        }
        else if (faceSrcFaces[fi] == -1 && faceTgtFaces[fi] != -1)
        {
            const label tgtFacei = tgtPatchStart + faceTgtFaces[fi];
            const label tgtOwn = mesh_.faceOwner()[tgtFacei];

            // Get target face zone info
            const label zoneIndex = mesh_.faceZones().whichZone(tgtFacei);
            bool zoneFlip = false;
            if (zoneIndex >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneIndex];
                zoneFlip = fZone.flipMap()[fZone.whichFace(tgtFacei)];
            }

            meshMod.addFace
            (
                f,          // Face to add
                tgtOwn,     // Owner cell
                -1,         // Neighbour cell
                -1,         // Master point index
                -1,         // Master edge index
                tgtFacei,   // Master face index
                false,      // Flip
                tgtPatchi,  // Patch index
                zoneIndex,     // Zone index
                zoneFlip    // Zone sign
            );

            if (debug)
            {
                Info<< "Adding patch face " << f
                    << " owner celli:" << tgtOwn
                    << " patchi:" << tgtPatchi << endl;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Both faceSrcFaces and faceTgtFaces are -1 for face " << fi
                << exit(FatalError);
        }
    }
}


void Foam::mergePatchPairs::addEdgeAddedPoints
(
    HashTable<labelList, edge, Hash<edge>>& edgeAddedPoints,
    const primitivePatch& patch,
    const List<DynamicList<label>>& patchEdgeAddedPoints
) const
{
    const edgeList& patchEdges = patch.edges();
    const labelList& pmp = patch.meshPoints();

    forAll(patchEdgeAddedPoints, ei)
    {
        if (patchEdgeAddedPoints[ei].size())
        {
            labelList newEdgePoints(patchEdgeAddedPoints[ei].size());

            forAll(newEdgePoints, npi)
            {
                newEdgePoints[npi] = newPoints_[patchEdgeAddedPoints[ei][npi]];
            }

            edgeAddedPoints.insert
            (
                edge(pmp[patchEdges[ei].start()], pmp[patchEdges[ei].end()]),
                newEdgePoints
            );
        }
    }
}


void Foam::mergePatchPairs::updatePoints
(
    labelList& pointMap,
    const primitivePatch& patch,
    const labelList& pointPoints
) const
{
    const labelList& pmp = patch.meshPoints();

    forAll(pointPoints, ppi)
    {
        pointMap[pmp[ppi]] = newPoints_[pointPoints[ppi]];
    }
}


void Foam::mergePatchPairs::modifyFaces
(
    polyTopoChange& meshMod,
    const polyPatchIntersection& intersection
) const
{
    HashTable<labelList, edge, Hash<edge>> edgeAddedPoints;
    addEdgeAddedPoints
    (
        edgeAddedPoints,
        intersection.srcPatch(),
        intersection.srcEdgePoints()
    );
    addEdgeAddedPoints
    (
        edgeAddedPoints,
        intersection.tgtPatch(),
        intersection.tgtEdgePoints()
    );

    labelList pointMap(identityMap(mesh_.nPoints()));
    updatePoints
    (
        pointMap,
        intersection.srcPatch(),
        intersection.srcPointPoints()
    );
    updatePoints
    (
        pointMap,
        intersection.tgtPatch(),
        intersection.tgtPointPoints()
    );

    DynamicList<label> modifiedFace;

    const faceList& meshFaces = mesh_.faces();

    forAll(meshFaces, fi)
    {
        if (meshMod.faceRemoved(fi))
        {
            continue;
        }

        const face& f = meshFaces[fi];

        bool addNext = true;
        bool modified = false;

        forAll(f, pi)
        {
            const edge e(f[pi], f[f.fcIndex(pi)]);

            const HashTable<labelList, edge, Hash<edge>>::const_iterator iter =
                edgeAddedPoints.find(e);

            if (iter != edgeAddedPoints.end())
            {
                const labelList& addedPoints = iter();

                if (edge::compare(e, iter.key()) == 1)
                {
                    modifiedFace.append(addedPoints);
                }
                else
                {
                    forAllReverse(addedPoints, i)
                    {
                        modifiedFace.append(addedPoints[i]);
                    }
                }

                if (f[f.fcIndex(pi)] == f[0])
                {
                    // modifiedFace.first() = modifiedFace.last();
                    // modifiedFace.setSize(modifiedFace.size() - 1);
                    modifiedFace.erase(modifiedFace.begin());
                }

                addNext = false;
                modified = true;
            }
            else
            {
                if (addNext)
                {
                    modifiedFace.append(pointMap[f[pi]]);
                }
                else
                {
                    addNext = true;
                }
            }
        }

        // Update points of point-connected faces if necessary
        if (!modified)
        {
            forAll(f, pi)
            {
                if (pointMap[f[pi]] != f[pi])
                {
                    if (!modified)
                    {
                        modifiedFace = f;
                        modified = true;
                    }

                    modifiedFace[pi] = pointMap[f[pi]];
                }
            }
        }

        if (modified)
        {
            // Get current zone info
            const label zoneIndex = mesh_.faceZones().whichZone(fi);
            bool zoneFlip = false;
            if (zoneIndex >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneIndex];
                zoneFlip = fZone.flipMap()[fZone.whichFace(fi)];
            }

            if (mesh_.isInternalFace(fi))
            {
                if (debug)
                {
                    Info<< "Modifying internal face " << fi
                        << " newface:" << modifiedFace
                        << " oldFace:" << f
                        << " centre:" << mesh_.faceCentres()[fi]
                        << " owner:" << mesh_.faceOwner()[fi]
                        << " neighbour:" << mesh_.faceNeighbour()[fi]
                        << endl;
                }

                meshMod.modifyFace
                (
                    face(modifiedFace),     // Modified face
                    fi,                     // index of face being modified
                    mesh_.faceOwner()[fi],  // Owner cell
                    mesh_.faceNeighbour()[fi], // Neighbour cell
                    false,                  // Face flip
                    -1,                     // Patch index
                    zoneIndex,                 // Zone index
                    zoneFlip                // Zone flip
                );
            }
            else
            {
                if (debug)
                {
                    Info<< "Modifying boundary face " << fi
                        << " newface:" << modifiedFace
                        << " oldFace:" << f
                        << " centre:" << mesh_.faceCentres()[fi]
                        << " owner:" << mesh_.faceOwner()[fi]
                        << " patch:" << mesh_.boundaryMesh().whichPatch(fi)
                        << endl;
                }

                meshMod.modifyFace
                (
                    face(modifiedFace),         // Modified face
                    fi,                         // index of face being modified
                    mesh_.faceOwner()[fi],      // Owner cell
                    -1,                         // Neighbour cell
                    false,                      // Face flip
                    mesh_.boundaryMesh().whichPatch(fi), // Patch index
                    zoneIndex,                     // Zone index
                    zoneFlip                    // Zone flip
                );
            }
        }

        modifiedFace.clear();
    }
}


void Foam::mergePatchPairs::intersectPatchPair
(
    polyTopoChange& meshMod,
    const polyPatch& srcPatch,
    const polyPatch& tgtPatch
) const
{
    polyPatchIntersection intersection(srcPatch, tgtPatch, snapTol_);
    removeFaces(meshMod, intersection);
    addPoints(meshMod, intersection);
    addFaces(meshMod, intersection);
    modifyFaces(meshMod, intersection);
}


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::mergePatchPairs::merge
(
    const List<Pair<label>>& patchPairs
) const
{
    // Construct the mesh modifier for merging patch-pairs
    polyTopoChange meshMod(mesh_);

    removePoints(meshMod, patchPairs);

    // Accumulate the mesh modifications for merging the patch-pairs
    forAll(patchPairs, ppi)
    {
        // Intersect patch-pair and update the mesh modifier
        intersectPatchPair
        (
            meshMod,
            mesh_.boundaryMesh()[patchPairs[ppi].first()],
            mesh_.boundaryMesh()[patchPairs[ppi].second()]
        );
    }

    // Change mesh and return map
    return meshMod.changeMesh
    (
        mesh_,
        false // Update mesh without moving points and calculating meshPhi
    );
}


bool Foam::mergePatchPairs::connected
(
    boolList& patchPoints,
    const label patchi
) const
{
    const labelList& pmp = mesh_.boundaryMesh()[patchi].meshPoints();

    forAll(pmp, pi)
    {
        if (patchPoints[pmp[pi]])
        {
            return true;
        }

        patchPoints[pmp[pi]] = true;
    }

    return false;
}


template<class Mesh>
inline void Foam::mergePatchPairs::merge
(
    Mesh& mesh,
    const List<Pair<word>>& patchPairNames
) const
{
    const List<Pair<label>> patchPairs(this->patchPairs(patchPairNames));

    bool connected = false;

    if (patchPairs.size() > 1)
    {
        boolList patchPoints(mesh_.nPoints(), false);

        forAll(patchPairs, ppi)
        {
            if
            (
                this->connected(patchPoints, patchPairs[ppi].first())
             || this->connected(patchPoints, patchPairs[ppi].second())
            )
            {
                connected = true;
                break;
            }
        }
    }

    if (connected)
    {
        // Merge patch-pairs and update fields
        forAll(patchPairNames, ppi)
        {
            Info<< "Merging patch pair " << patchPairNames[ppi] << endl;

            mesh.topoChange
            (
                merge(List<Pair<label>>(patchPairs[ppi]))
            );
        }
    }
    else
    {
        Info<< "Merging patch pairs " << patchPairNames << endl;

        mesh.topoChange(merge(patchPairs));
    }

    // Find the stitched patches that are now zero-sized
    labelHashSet zeroSizedMergedPatches;

    forAll(patchPairs, ppi)
    {
        if (mesh.boundaryMesh()[patchPairs[ppi].first()].size() == 0)
        {
            zeroSizedMergedPatches.insert(patchPairs[ppi].first());
        }

        if (mesh.boundaryMesh()[patchPairs[ppi].second()].size() == 0)
        {
            zeroSizedMergedPatches.insert(patchPairs[ppi].second());
        }
    }

    // Create a list of the remaining patch old patch labels
    labelList remainingPatches
    (
        mesh.boundaryMesh().size() - zeroSizedMergedPatches.size()
    );

    label rpi = 0;
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (!zeroSizedMergedPatches.found(patchi))
        {
            remainingPatches[rpi++] = patchi;
        }
    }

    // Remove the zero-sized stitched patches
    mesh_.reorderPatches(remainingPatches, true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mergePatchPairs::mergePatchPairs
(
    polyMesh& mesh,
    const List<Pair<word>>& patchPairNames,
    const scalar snapTol
)
:
    mesh_(mesh),
    snapTol_(snapTol)
{
    // Merge patch-pairs and update fields
    merge(mesh, patchPairNames);
}


Foam::mergePatchPairs::mergePatchPairs
(
    fvMesh& mesh,
    const List<Pair<word>>& patchPairNames,
    const scalar snapTol
)
:
    mesh_(mesh),
    snapTol_(snapTol)
{
    // Merge patch-pairs and update fields
    merge(mesh, patchPairNames);
}


// ************************************************************************* //
