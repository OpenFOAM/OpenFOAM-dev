/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "nonConformalBoundary.H"
#include "nonConformalCoupledPolyPatch.H"
#include "nonConformalErrorPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class FaceList, class PointField>
edge meshEdge
(
    const PrimitivePatch<FaceList, PointField>& p,
    const label edgei
)
{
    return edge
    (
        p.meshPoints()[p.edges()[edgei][0]],
        p.meshPoints()[p.edges()[edgei][1]]
    );
}

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalBoundary, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::indirectPrimitivePatch
Foam::nonConformalBoundary::boundary(const labelList& patches) const
{
    DynamicList<label> faces;
    forAll(patches, i)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patches[i]];

        faces.append(identity(pp.size()) + pp.start());
    }

    return
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh().faces(), faces),
            mesh().points()
        );
}


const Foam::labelList&
Foam::nonConformalBoundary::meshPointOwnerOrigBoundaryPoint() const
{
    if (!meshPointOwnerOrigBoundaryPointPtr_.valid())
    {
        meshPointOwnerOrigBoundaryPointPtr_.set
        (
            new labelList(mesh().nPoints(), -1)
        );

        forAll(ownerOrigBoundary_.meshPoints(), ownerOrigBoundaryPointi)
        {
            meshPointOwnerOrigBoundaryPointPtr_()
                [ownerOrigBoundary_.meshPoints()[ownerOrigBoundaryPointi]] =
                ownerOrigBoundaryPointi;
        }
    }

    return meshPointOwnerOrigBoundaryPointPtr_();
}


const Foam::vectorField&
Foam::nonConformalBoundary::ownerOrigBoundaryPointNormals() const
{
    if (!ownerOrigBoundaryPointNormalsPtr_.valid())
    {
        const faceList& faces = ownerOrigBoundary_.localFaces();
        const vectorField faceNormals = ownerOrigBoundary_.faceNormals();

        vectorField pointNormals(ownerOrigBoundary_.nPoints(), Zero);

        forAll(faces, facei)
        {
            forAll(faces[facei], facePointi)
            {
                pointNormals[faces[facei][facePointi]] += faceNormals[facei];
            }
        }

        syncTools::syncPointList
        (
            mesh(),
            ownerOrigBoundary_.meshPoints(),
            pointNormals,
            plusEqOp<vector>(),
            vector::zero
        );

        ownerOrigBoundaryPointNormalsPtr_.set
        (
            (pointNormals/(mag(pointNormals) + vSmall)).ptr()
        );
    }

    return ownerOrigBoundaryPointNormalsPtr_();
}


const Foam::vectorField&
Foam::nonConformalBoundary::ownerOrigBoundaryPointNormals0() const
{
    if (!ownerOrigBoundaryPointNormals0Ptr_.valid())
    {
        const faceList& faces = ownerOrigBoundary_.localFaces();

        vectorField pointNormals(ownerOrigBoundary_.nPoints(), Zero);

        forAll(faces, facei)
        {
            forAll(faces[facei], facePointi)
            {
                pointNormals[faces[facei][facePointi]] +=
                    faces[facei].normal(mesh().oldPoints());
            }
        }

        syncTools::syncPointList
        (
            mesh(),
            ownerOrigBoundary_.meshPoints(),
            pointNormals,
            plusEqOp<vector>(),
            vector::zero
        );

        ownerOrigBoundaryPointNormals0Ptr_.set
        (
            (pointNormals/(mag(pointNormals) + vSmall)).ptr()
        );
    }

    return ownerOrigBoundaryPointNormals0Ptr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalBoundary::nonConformalBoundary(const polyMesh& mesh)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        nonConformalBoundary
    >(mesh),
    ownerOrigBoundary_(boundary(ownerOrigPatchIDs())),
    meshPointOwnerOrigBoundaryPointPtr_(nullptr),
    ownerOrigBoundaryPointMeshPointPtr_(nullptr),
    ownerOrigBoundaryEdgeMeshEdgePtr_(nullptr),
    ownerOrigBoundaryEdgesPtr_(nullptr),
    ownerOrigBoundaryMeshEdgesPtr_(nullptr),
    patchPointOwnerOrigBoundaryPointsPtr_(mesh.boundaryMesh().size()),
    patchEdgeOwnerOrigBoundaryEdgesPtr_(mesh.boundaryMesh().size()),
    ownerOrigBoundaryPointNormalsPtr_(nullptr)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalBoundary::~nonConformalBoundary()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::nonConformalBoundary::movePoints()
{
    ownerOrigBoundary_.clearGeom();

    ownerOrigBoundaryPointNormalsPtr_.clear();

    return true;
}


Foam::labelList Foam::nonConformalBoundary::nonConformalNonCoupledPatchIDs
(
    const label side,
    label (nonConformalCoupledPolyPatch::*method)() const
) const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    labelHashSet origPatchIDTable;
    DynamicList<label> nonCoupledPatchIDs(pbm.size());

    forAll(pbm, nccPatchi)
    {
        const polyPatch& pp = pbm[nccPatchi];

        if (isA<nonConformalCoupledPolyPatch>(pp))
        {
            const nonConformalCoupledPolyPatch& nccPp =
                refCast<const nonConformalCoupledPolyPatch>(pbm[nccPatchi]);

            if (side == 1 && !nccPp.owner()) continue;

            if (side == -1 && !nccPp.neighbour()) continue;

            if (origPatchIDTable.found(nccPp.origPatchID())) continue;

            origPatchIDTable.insert(nccPp.origPatchID());
            nonCoupledPatchIDs.append((nccPp.*method)());
        }
    }

    labelList result;
    result.transfer(nonCoupledPatchIDs);
    return result;
}


Foam::labelList Foam::nonConformalBoundary::allOrigPatchIDs() const
{
    auto method = &nonConformalCoupledPolyPatch::origPatchID;
    return nonConformalNonCoupledPatchIDs(0, method);
}


Foam::labelList Foam::nonConformalBoundary::allErrorPatchIDs() const
{
    auto method = &nonConformalCoupledPolyPatch::errorPatchID;
    return nonConformalNonCoupledPatchIDs(0, method);
}


Foam::labelList Foam::nonConformalBoundary::ownerOrigPatchIDs() const
{
    auto method = &nonConformalCoupledPolyPatch::origPatchID;
    return nonConformalNonCoupledPatchIDs(1, method);
}


Foam::labelList Foam::nonConformalBoundary::ownerErrorPatchIDs() const
{
    auto method = &nonConformalCoupledPolyPatch::errorPatchID;
    return nonConformalNonCoupledPatchIDs(1, method);
}


const Foam::labelList&
Foam::nonConformalBoundary::ownerOrigBoundaryPointMeshPoint() const
{
    if (!ownerOrigBoundaryPointMeshPointPtr_.valid())
    {
        ownerOrigBoundaryPointMeshPointPtr_.set
        (
            new labelList(ownerOrigBoundary_.meshPoints())
        );

        // ...
        meshPointOwnerOrigBoundaryPoint();
        labelList& map = meshPointOwnerOrigBoundaryPointPtr_();

        // ...
        label ownerOrigBoundaryPointi = ownerOrigBoundary_.nPoints();
        DynamicList<label> remotePoints;

        // ...
        for
        (
            label ownerOrigBoundaryEdgei = ownerOrigBoundary_.nEdges();
            ownerOrigBoundaryEdgei < ownerOrigBoundaryEdgeMeshEdge().size();
            ++ ownerOrigBoundaryEdgei
        )
        {
            const label meshEdgei =
                ownerOrigBoundaryEdgeMeshEdge()[ownerOrigBoundaryEdgei];

            const edge& e = mesh().edges()[meshEdgei];

            forAll(e, i)
            {
                const label meshPointi = e[i];

                if (map[meshPointi] == -1)
                {
                    map[meshPointi] = ownerOrigBoundaryPointi ++;
                    remotePoints.append(meshPointi);
                }
            }
        }

        // ...
        ownerOrigBoundaryPointMeshPointPtr_->append(remotePoints);
    }

    return ownerOrigBoundaryPointMeshPointPtr_();
}


const Foam::labelList&
Foam::nonConformalBoundary::ownerOrigBoundaryEdgeMeshEdge() const
{
    if (!ownerOrigBoundaryEdgeMeshEdgePtr_.valid())
    {
        // Create boundary of all owner-orig and proc patches
        labelList ownerOrigAndProcPatchIDs = this->ownerOrigPatchIDs();
        forAll(mesh().boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(mesh().boundaryMesh()[patchi]))
            {
                ownerOrigAndProcPatchIDs.append(patchi);
            }
        }
        const indirectPrimitivePatch ownerOrigAndProcBoundary
        (
            boundary(ownerOrigAndProcPatchIDs)
        );

        // Create the mesh edge mapping for the owner-orig-and-proc boundary
        labelList ownerOrigAndProcBoundaryMeshEdges
        (
            ownerOrigAndProcBoundary.meshEdges
            (
                mesh().edges(),
                mesh().pointEdges()
            )
        );

        // Map from owner-orig to owner-orig-and-proc edge
        const labelList map =
            primitivePatch
            (
                SubList<face>
                (
                    ownerOrigAndProcBoundary.localFaces(),
                    ownerOrigBoundary_.size()
                ),
                ownerOrigAndProcBoundary.localPoints()
            )
           .meshEdges
            (
                ownerOrigAndProcBoundary.edges(),
                ownerOrigAndProcBoundary.pointEdges()
            );

        // Initialise the mapping for the owner-orig boundary. Note, this only
        // connects edges that are part of the local owner-orig boundary. Edges
        // that are connected to the owner-orig boundary on a different process
        // are added below.
        ownerOrigBoundaryEdgeMeshEdgePtr_.set
        (
            new labelList
            (
                labelField(ownerOrigAndProcBoundaryMeshEdges, map)
            )
        );

        // Create a reverse map from owner-orig-and-proc to owner-orig edge
        labelList rMap(ownerOrigAndProcBoundary.nEdges(), -1);
        forAll(map, i)
        {
            rMap[map[i]] = i;
        }

        // Synchronise
        syncTools::syncEdgeList
        (
            mesh(),
            ownerOrigAndProcBoundaryMeshEdges,
            rMap,
            maxEqOp<label>(),
            label(-1)
        );

        // Remove all local indexing from the map
        forAll(map, i)
        {
            rMap[map[i]] = -1;
        }

        // Any valid rMap entries now refer to edges that are remotely
        // connected to the owner-orig boundary. Append them to the addressing.
        DynamicList<label> remoteMeshEdges;
        forAll(rMap, i)
        {
            if (rMap[i] != -1)
            {
                remoteMeshEdges.append(ownerOrigAndProcBoundaryMeshEdges[i]);
            }
        }

        ownerOrigBoundaryEdgeMeshEdgePtr_->append(remoteMeshEdges);
    }

    return ownerOrigBoundaryEdgeMeshEdgePtr_();
}


const Foam::edgeList&
Foam::nonConformalBoundary::ownerOrigBoundaryEdges() const
{
    if (!ownerOrigBoundaryEdgesPtr_.valid())
    {
        ownerOrigBoundaryEdgesPtr_.set
        (
            new edgeList(ownerOrigBoundary_.edges())
        );

        const labelList& map = meshPointOwnerOrigBoundaryPoint();

        DynamicList<edge> remoteEdges;

        for
        (
            label ownerOrigBoundaryEdgei = ownerOrigBoundary_.nEdges();
            ownerOrigBoundaryEdgei < ownerOrigBoundaryEdgeMeshEdge().size();
            ++ ownerOrigBoundaryEdgei
        )
        {
            const label meshEdgei =
                ownerOrigBoundaryEdgeMeshEdge()[ownerOrigBoundaryEdgei];

            const edge& e = mesh().edges()[meshEdgei];

            remoteEdges.append(edge(map[e.start()], map[e.end()]));
        }

        ownerOrigBoundaryEdgesPtr_->append(remoteEdges);
    }

    return ownerOrigBoundaryEdgesPtr_();
}


const Foam::edgeList&
Foam::nonConformalBoundary::ownerOrigBoundaryMeshEdges() const
{
    if (!ownerOrigBoundaryMeshEdgesPtr_.valid())
    {
        const edgeList& edges = ownerOrigBoundaryEdges();

        const labelList& pointMeshPoint =
            ownerOrigBoundaryPointMeshPoint();

        ownerOrigBoundaryMeshEdgesPtr_.set
        (
            new edgeList(edges.size())
        );

        forAll(edges, edgei)
        {
            ownerOrigBoundaryMeshEdgesPtr_()[edgei] =
                edge
                (
                    pointMeshPoint[edges[edgei].start()],
                    pointMeshPoint[edges[edgei].end()]
                );
        }
    }

    return ownerOrigBoundaryMeshEdgesPtr_();
}


const Foam::labelList&
Foam::nonConformalBoundary::patchPointOwnerOrigBoundaryPoints
(
    const label patchi
) const
{
    if (!patchPointOwnerOrigBoundaryPointsPtr_.set(patchi))
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        const faceList patchOwnerOrigBoundaryLocalFaces
        (
            renumber
            (
                meshPointOwnerOrigBoundaryPoint(),
                static_cast<const faceList&>(pp)
            )
        );

        const primitivePatch ownerOrigBoundaryLocalPatch
        (
            SubList<face>(patchOwnerOrigBoundaryLocalFaces, pp.size()),
            ownerOrigBoundary_.localPoints()
        );

        patchPointOwnerOrigBoundaryPointsPtr_.set
        (
            patchi,
            new labelList(ownerOrigBoundaryLocalPatch.meshPoints())
        );

        patchEdgeOwnerOrigBoundaryEdgesPtr_.set
        (
            patchi,
            new labelList
            (
                ownerOrigBoundaryLocalPatch.meshEdges
                (
                    ownerOrigBoundary_.edges(),
                    ownerOrigBoundary_.pointEdges()
                )
            )
        );

        // Check the addressing is valid
        if (debug)
        {
            const labelList& ppPointOwnerOrigBoundaryPoints =
                patchPointOwnerOrigBoundaryPointsPtr_[patchi];

            forAll(pp.meshPoints(), ppPointi)
            {
                const label ownerOrigBoundaryPointi =
                    ppPointOwnerOrigBoundaryPoints[ppPointi];

                if
                (
                    pp.meshPoints()[ppPointi]
                 != ownerOrigBoundary_.meshPoints()[ownerOrigBoundaryPointi]
                )
                {
                    FatalErrorInFunction
                        << "Patch point does not match all boundary point"
                        << exit(FatalError);
                }
            }

            const labelList& ppEdgeOwnerOrigBoundaryEdges =
                patchEdgeOwnerOrigBoundaryEdgesPtr_[patchi];

            forAll(pp.edges(), ppEdgei)
            {
                const label ownerOrigBoundaryEdgei =
                    ppEdgeOwnerOrigBoundaryEdges[ppEdgei];

                if
                (
                    meshEdge(pp, ppEdgei)
                 != meshEdge(ownerOrigBoundary_, ownerOrigBoundaryEdgei)
                )
                {
                    FatalErrorInFunction
                        << "Patch edge does not match all boundary edge"
                        << exit(FatalError);
                }
            }
        }
    }

    return patchPointOwnerOrigBoundaryPointsPtr_[patchi];
}


const Foam::labelList&
Foam::nonConformalBoundary::patchEdgeOwnerOrigBoundaryEdges
(
    const label patchi
) const
{
    if (!patchEdgeOwnerOrigBoundaryEdgesPtr_.set(patchi))
    {
        patchPointOwnerOrigBoundaryPoints(patchi);
    }

    return patchEdgeOwnerOrigBoundaryEdgesPtr_[patchi];
}


Foam::tmp<Foam::vectorField>
Foam::nonConformalBoundary::patchPointNormals(const label patchi) const
{
    return
        tmp<vectorField>
        (
            new vectorField
            (
                UIndirectList<vector>
                (
                    ownerOrigBoundaryPointNormals(),
                    patchPointOwnerOrigBoundaryPoints(patchi)
                )
            )
        );
}


Foam::tmp<Foam::vectorField>
Foam::nonConformalBoundary::patchPointNormals0(const label patchi) const
{
    return
        tmp<vectorField>
        (
            new vectorField
            (
                UIndirectList<vector>
                (
                    ownerOrigBoundaryPointNormals0(),
                    patchPointOwnerOrigBoundaryPoints(patchi)
                )
            )
        );
}


// ************************************************************************* //
