/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshStitcher.H"
#include "globalIndex.H"
#include "fvcSurfaceIntegrate.H"
#include "meshObjects.H"
#include "polyTopoChangeMap.H"
#include "syncTools.H"
#include "surfaceToVolVelocity.H"

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
    defineTypeNameAndDebug(fvMeshStitcher, 0);
    defineRunTimeSelectionTable(fvMeshStitcher, fvMesh);
}


const Foam::word Foam::fvMeshStitcher::nccFieldPrefix_ =
    fvMeshStitcher::typeName + ":";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshStitcher::intersectNonConformalCyclic
(
    const nonConformalCyclicFvPatch& nccFvp,
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf,
    const tmp<surfaceLabelField::Boundary>& tOrigFacesNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigSfNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigCfNbrBf,
    List<part>& origEdgeParts
) const
{
    const nonConformalCyclicFvPatch& nbrNccFvp = nccFvp.nbrPatch();

    // Alias the original poly patches
    const polyPatch& origPp = nccFvp.origPatch().patch();
    const polyPatch& nbrOrigPp = nbrNccFvp.origPatch().patch();

    // Get the indices of the non-conformal patches
    labelList patchis(Pstream::nProcs(), -1);
    labelList nbrPatchis(Pstream::nProcs(), -1);
    patchis[Pstream::myProcNo()] = nccFvp.index();
    nbrPatchis[Pstream::myProcNo()] = nbrNccFvp.index();
    forAll(mesh_.boundary(), patchj)
    {
        const fvPatch& fvp = mesh_.boundary()[patchj];

        if (isA<nonConformalProcessorCyclicFvPatch>(fvp))
        {
            const nonConformalProcessorCyclicFvPatch& ncpcFvp =
                refCast<const nonConformalProcessorCyclicFvPatch>(fvp);

            if (ncpcFvp.referPatchID() == nccFvp.index())
            {
                patchis[ncpcFvp.neighbProcNo()] = patchj;
            }
            if (ncpcFvp.referPatchID() == nbrNccFvp.index())
            {
                nbrPatchis[ncpcFvp.neighbProcNo()] = patchj;
            }
        }
    }

    // Get the intersection geometry
    const patchToPatches::intersection& intersection =
        nccFvp.nonConformalCyclicPatch().intersection();

    // Unpack the patchToPatch addressing into lists of indices (fixed lists of
    // 3 labels; owner face, neighbour face, couple index). These will be used
    // to create the non-conformal faces, so sort them to make sure the
    // non-conformal interfaces are ordered.
    auto procFacesToIndices = []
    (
        const List<List<remote>>& faceOtherProcFaces,
        const bool owner
    )
    {
        // Compute interface sizes
        labelList is(Pstream::nProcs(), 0);
        forAll(faceOtherProcFaces, facei)
        {
            forAll(faceOtherProcFaces[facei], i)
            {
                const remote& otherProcFacei = faceOtherProcFaces[facei][i];

                is[otherProcFacei.proci] ++;
            }
        }

        // Allocate the indices
        List<List<FixedList<label, 3>>> indices(Pstream::nProcs());
        forAll(indices, proci)
        {
            indices[proci].resize(is[proci]);
        }

        // Set the indices
        is = 0;
        forAll(faceOtherProcFaces, facei)
        {
            forAll(faceOtherProcFaces[facei], i)
            {
                const remote& otherProcFacei = faceOtherProcFaces[facei][i];

                FixedList<label, 3>& index =
                    indices[otherProcFacei.proci][is[otherProcFacei.proci] ++];

                index = {facei, otherProcFacei.elementi, i};

                if (!owner) Swap(index[0], index[1]);
            }
        }

        // Sort to ensure ordering
        forAll(indices, proci)
        {
            sort(indices[proci]);
        }

        return indices;
    };

    List<List<FixedList<label, 3>>> indices =
        procFacesToIndices(intersection.srcTgtProcFaces(), true);

    List<List<FixedList<label, 3>>> nbrIndices =
        procFacesToIndices(intersection.tgtSrcProcFaces(), false);

    // If addressing has been provided, then modify the indices to match. When
    // a coupling has to be added, the couple index is set to -1. This
    // indicates that there is no geometry in the patchToPatch engine for this
    // coupling, and for a small stabilisation value to be used instead.
    auto matchIndices = [&polyFacesBf, &tOrigFacesNbrBf]
    (
        List<List<FixedList<label, 3>>>& indices,
        const nonConformalCyclicFvPatch& nccFvp,
        const polyPatch& origPp,
        const labelList& patchis,
        const bool owner
    )
    {
        const surfaceLabelField::Boundary& origFacesNbrBf = tOrigFacesNbrBf();

        // Create what the indices should be
        List<List<FixedList<label, 3>>> indicesRef(indices.size());
        forAll(indices, proci)
        {
            const label patchi = patchis[proci];

            if (patchi != -1)
            {
                indicesRef[proci].resize(polyFacesBf[patchi].size());

                forAll(polyFacesBf[patchi], patchFacei)
                {
                    FixedList<label, 3>& indexRef =
                        indicesRef[proci][patchFacei];

                    indexRef =
                        {
                            polyFacesBf[patchi][patchFacei] - origPp.start(),
                            origFacesNbrBf[patchi][patchFacei],
                            -1
                        };

                    if (!owner) Swap(indexRef[0], indexRef[1]);
                }
            }
        }

        // Populate the indices with the index of the coupling
        label nCouplesRemoved = 0, nCouplesAdded = 0;
        forAll(indices, proci)
        {
            const label patchi = patchis[proci];

            if (patchi != -1)
            {
                label patchFacei = 0, i = 0;

                while
                (
                    patchFacei < polyFacesBf[patchi].size()
                 && i < indices[proci].size()
                )
                {
                    const FixedList<label, 3> index
                    ({
                        indices[proci][i][0],
                        indices[proci][i][1],
                        -1
                    });

                    FixedList<label, 3>& indexRef =
                        indicesRef[proci][patchFacei];

                    if (index < indexRef)
                    {
                        nCouplesRemoved ++;
                        i ++;
                    }
                    else if (index == indexRef)
                    {
                        indexRef[2] = indices[proci][i][2];
                        patchFacei ++;
                        i ++;
                    }
                    else // (index > indexRef)
                    {
                        nCouplesAdded ++;
                        patchFacei ++;
                    }
                }

                nCouplesRemoved +=
                    min(indices[proci].size() - i, 0);
                nCouplesAdded +=
                    min(polyFacesBf[patchi].size() - patchFacei, 0);
            }
        }

        // Report if changes have been made
        reduce(nCouplesRemoved, sumOp<label>());
        reduce(nCouplesAdded, sumOp<label>());
        if ((nCouplesRemoved || nCouplesAdded) && nccFvp.owner())
        {
            Info<< indent << nCouplesRemoved << '/' << nCouplesAdded
                << " small couplings removed/added to " << nccFvp.name()
                << endl;
        }

        // Set the indices to the correct values
        Swap(indices, indicesRef);
    };

    if (tOrigFacesNbrBf.valid())
    {
        matchIndices(indices, nccFvp, origPp, patchis, true);
        matchIndices(nbrIndices, nbrNccFvp, nbrOrigPp, nbrPatchis, false);
    }

    // Create couplings by transferring geometry from the original to the
    // non-conformal patches
    auto createCouplings =
    [
        &polyFacesBf,
        &tOrigSfNbrBf,
        &tOrigCfNbrBf,
        &SfBf,
        &CfBf
    ]
    (
        const List<List<FixedList<label, 3>>>& indices,
        const List<DynamicList<couple>>& couples,
        const nonConformalCyclicFvPatch& nccFvp,
        const polyPatch& origPp,
        const labelList& patchis,
        const bool owner
    )
    {
        forAll(patchis, proci)
        {
            const label patchi = patchis[proci];
            const label patchSize = indices[proci].size();

            if (patchi == -1 && patchSize)
            {
                FatalErrorInFunction
                    << "Intersection generated coupling through "
                    << nccFvp.name() << " between processes "
                    << Pstream::myProcNo() << " and " << proci
                    << " but a corresponding processor interface was not found"
                    << exit(FatalError);
            }

            if (patchi != -1)
            {
                polyFacesBf[patchi].resize(patchSize);
                SfBf[patchi].resize(patchSize);
                CfBf[patchi].resize(patchSize);

                forAll(indices[proci], patchFacei)
                {
                    const label origFacei = indices[proci][patchFacei][!owner];
                    const label i = indices[proci][patchFacei][2];

                    polyFacesBf[patchi][patchFacei] =
                        origFacei + origPp.start();

                    couple c;
                    if (i != -1)
                    {
                        c = couples[origFacei][i];
                    }
                    else
                    {
                        c =
                            couple
                            (
                                part
                                (
                                    small*origPp.faceAreas()[origFacei],
                                    origPp.faceCentres()[origFacei]
                                ),
                                part
                                (
                                    small*tOrigSfNbrBf()[patchi][patchFacei],
                                    tOrigCfNbrBf()[patchi][patchFacei]
                                )
                            );
                    }
                    if (!owner)
                    {
                        c.nbr = c;
                    }

                    SfBf[patchi][patchFacei] = c.nbr.area;
                    CfBf[patchi][patchFacei] = c.nbr.centre;

                    part origP
                    (
                        SfBf[origPp.index()][origFacei],
                        CfBf[origPp.index()][origFacei]
                    );
                    origP -= c;

                    SfBf[origPp.index()][origFacei] = origP.area;
                    CfBf[origPp.index()][origFacei] = origP.centre;
                }
            }
        }
    };

    createCouplings
    (
        indices,
        intersection.srcCouples(),
        nccFvp,
        origPp,
        patchis,
        true
    );

    createCouplings
    (
        nbrIndices,
        intersection.tgtCouples(),
        nbrNccFvp,
        nbrOrigPp,
        nbrPatchis,
        false
    );

    // Add error geometry to the original patches
    forAll(intersection.srcErrorParts(), origFacei)
    {
        part origP
        (
            SfBf[origPp.index()][origFacei],
            CfBf[origPp.index()][origFacei]
        );
        origP += intersection.srcErrorParts()[origFacei];

        SfBf[origPp.index()][origFacei] = origP.area;
        CfBf[origPp.index()][origFacei] = origP.centre;
    }

    // Store the orig edge parts
    forAll(intersection.srcEdgeParts(), origEdgei)
    {
        origEdgeParts[origEdgei] += intersection.srcEdgeParts()[origEdgei];
    }
}


Foam::List<Foam::patchToPatches::intersection::part>
Foam::fvMeshStitcher::calculateOwnerOrigBoundaryEdgeParts
(
    const List<List<part>>& patchEdgeParts
) const
{
    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList& ownerOrigBoundaryPointMeshPoint =
        ncb.ownerOrigBoundaryPointMeshPoint();
    const labelList& ownerOrigBoundaryEdgeMeshEdge =
        ncb.ownerOrigBoundaryEdgeMeshEdge();
    const edgeList& ownerOrigBoundaryEdges = ncb.ownerOrigBoundaryEdges();
    const edgeList& ownerOrigBoundaryMeshEdges =
        ncb.ownerOrigBoundaryMeshEdges();

    // Sum up boundary edge parts, being careful with signs
    labelList ownerOrigBoundaryEdgeNParts(ownerOrigBoundaryEdges.size(), 0);
    List<part> ownerOrigBoundaryEdgeParts(ownerOrigBoundaryEdges.size());
    forAll(patchEdgeParts, patchi)
    {
        if (patchEdgeParts[patchi].empty()) continue;

        const polyPatch& patch = pbMesh[patchi];

        const labelList& patchEdgeOwnerOrigBoundaryEdges =
            ncb.patchEdgeOwnerOrigBoundaryEdges(patchi);

        forAll(patchEdgeParts[patchi], patchEdgei)
        {
            const label ownerOrigBoundaryEdgei =
                patchEdgeOwnerOrigBoundaryEdges[patchEdgei];

            const label sign =
                edge::compare
                (
                    meshEdge(patch, patchEdgei),
                    ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                );

            const part& pU = patchEdgeParts[patchi][patchEdgei];
            const part p = sign > 0 ? pU : -pU;

            ownerOrigBoundaryEdgeNParts[ownerOrigBoundaryEdgei] ++;
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei] += p;
        }
    }

    // Give every point in the all boundary a global unique index
    const globalIndex globalOwnerOrigBoundaryPoints
    (
        ownerOrigBoundaryPointMeshPoint.size()
    );
    labelList ownerOrigBoundaryPointIndices
    (
        ownerOrigBoundaryPointMeshPoint.size()
    );
    forAll(ownerOrigBoundaryPointIndices, ownerOrigBoundaryPointi)
    {
        ownerOrigBoundaryPointIndices[ownerOrigBoundaryPointi] =
            globalOwnerOrigBoundaryPoints.toGlobal(ownerOrigBoundaryPointi);
    }
    syncTools::syncPointList
    (
        mesh_,
        ownerOrigBoundaryPointMeshPoint,
        ownerOrigBoundaryPointIndices,
        minEqOp<label>(),
        labelMax
    );

    // Synchronise the sign of the edge parts so that they can be
    // meaningfully summed. This is done using the sign of the global point
    // indices; if they are not in order, then reverse the part area.
    forAll(ownerOrigBoundaryEdgeParts, ownerOrigBoundaryEdgei)
    {
        const edge& e = ownerOrigBoundaryEdges[ownerOrigBoundaryEdgei];

        if
        (
            ownerOrigBoundaryPointIndices[e.start()]
          > ownerOrigBoundaryPointIndices[e.end()])
        {
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei] =
                - ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];
        }
    }

    // Average the boundary edge parts across couplings
    syncTools::syncEdgeList
    (
        mesh_,
        ownerOrigBoundaryEdgeMeshEdge,
        ownerOrigBoundaryEdgeNParts,
        plusEqOp<label>(),
        label(0)
    );
    syncTools::syncEdgeList
    (
        mesh_,
        ownerOrigBoundaryEdgeMeshEdge,
        ownerOrigBoundaryEdgeParts,
        plusEqOp<part>(),
        part(),
        []
        (
            const transformer& vt,
            const bool forward,
            List<part>& fld
        )
        {
            if (forward)
            {
                forAll(fld, i)
                {
                    if (magSqr(fld[i].area) != 0)
                    {
                        fld[i].area = vt.transform(fld[i].area);
                        fld[i].centre = vt.transformPosition(fld[i].centre);
                    }
                }
            }
            else
            {
                forAll(fld, i)
                {
                    if (magSqr(fld[i].area) != 0)
                    {
                        fld[i].area = vt.invTransform(fld[i].area);
                        fld[i].centre = vt.invTransformPosition(fld[i].centre);
                    }
                }
            }
        }
    );
    forAll(ownerOrigBoundaryEdgeParts, ownerOrigBoundaryEdgei)
    {
        if (ownerOrigBoundaryEdgeNParts[ownerOrigBoundaryEdgei] != 0)
        {
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei].area /=
                ownerOrigBoundaryEdgeNParts[ownerOrigBoundaryEdgei];
        }
    }

    // Put the edge parts back to their original orientations
    forAll(ownerOrigBoundaryEdgeParts, ownerOrigBoundaryEdgei)
    {
        const edge& e = ownerOrigBoundaryEdges[ownerOrigBoundaryEdgei];

        if
        (
            ownerOrigBoundaryPointIndices[e.start()]
          > ownerOrigBoundaryPointIndices[e.end()]
        )
        {
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei] =
                - ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];
        }
    }

    return ownerOrigBoundaryEdgeParts;
}



void Foam::fvMeshStitcher::applyOwnerOrigBoundaryEdgeParts
(
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf,
    const List<part>& ownerOrigBoundaryEdgeParts
) const
{
    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList ownerOrigPatchIDs = ncb.ownerOrigPatchIDs();
    const labelList& ownerOrigBoundaryEdgeMeshEdge =
        ncb.ownerOrigBoundaryEdgeMeshEdge();
    const edgeList& ownerOrigBoundaryMeshEdges =
        ncb.ownerOrigBoundaryMeshEdges();

    boolList patchIsOwnerOrig(pbMesh.size(), false);
    UIndirectList<bool>(patchIsOwnerOrig, ownerOrigPatchIDs) = true;

    // Count the number of original faces attached to each boundary edge
    labelList ownerOrigBoundaryEdgeNOrigFaces
    (
        ownerOrigBoundaryEdgeMeshEdge.size(),
        0
    );
    forAll(ownerOrigBoundaryEdgeMeshEdge, ownerOrigBoundaryEdgei)
    {
        const label meshEdgei =
            ownerOrigBoundaryEdgeMeshEdge[ownerOrigBoundaryEdgei];

        forAll(mesh_.edgeFaces()[meshEdgei], edgeFacei)
        {
            const label facei = mesh_.edgeFaces()[meshEdgei][edgeFacei];

            const label patchi =
                mesh_.isInternalFace(facei)
              ? -1 : pbMesh.patchID()[facei - mesh_.nInternalFaces()];

            if (patchi != -1 && patchIsOwnerOrig[patchi])
            {
                ownerOrigBoundaryEdgeNOrigFaces[ownerOrigBoundaryEdgei] ++;
            }
        }
    }

    // Synchronise the boundary edge original face count
    syncTools::syncEdgeList
    (
        mesh_,
        ownerOrigBoundaryEdgeMeshEdge,
        ownerOrigBoundaryEdgeNOrigFaces,
        plusEqOp<label>(),
        label(0)
    );

    // If debugging, check that face changes are in sync
    tmp<surfaceLabelField> tChanged =
        debug
      ? surfaceLabelField::New("changed", mesh_, dimensioned<label>(dimless, 0))
      : tmp<surfaceLabelField>(nullptr);

    // Modify faces connected to the boundary edges
    forAll(ownerOrigBoundaryEdgeMeshEdge, ownerOrigBoundaryEdgei)
    {
        const label meshEdgei =
            ownerOrigBoundaryEdgeMeshEdge[ownerOrigBoundaryEdgei];

        const part& ownerOrigBoundaryEdgeP =
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];

        switch (ownerOrigBoundaryEdgeNOrigFaces[ownerOrigBoundaryEdgei])
        {
            case 0:
            {
                continue;
            }

            // If this edge is connected to just one original face then add any
            // edge area into that face
            case 1:
            {
                label origFacei = -1, origPatchi = -1, origPatchFacei = -1;
                forAll(mesh_.edgeFaces()[meshEdgei], edgeFacei)
                {
                    const label facei = mesh_.edgeFaces()[meshEdgei][edgeFacei];

                    const label patchi =
                        mesh_.isInternalFace(facei)
                      ? -1 : pbMesh.patchID()[facei - mesh_.nInternalFaces()];

                    if (patchi != -1 && patchIsOwnerOrig[patchi])
                    {
                        origFacei = facei;
                        origPatchi = patchi;
                        origPatchFacei = facei - pbMesh[patchi].start();
                        break;
                    }
                }

                if (origFacei != -1)
                {
                    vector& area =
                        SfSf.boundaryFieldRef()[origPatchi][origPatchFacei];
                    point& centre =
                        CfSf.boundaryFieldRef()[origPatchi][origPatchFacei];

                    const label sign =
                        mesh_.faces()[origFacei].edgeDirection
                        (
                            ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                        );

                    part p(area, centre);
                    p +=
                        sign > 0
                      ? ownerOrigBoundaryEdgeP
                      : -ownerOrigBoundaryEdgeP;

                    area = p.area;
                    centre = p.centre;

                    if (debug)
                    {
                        tChanged->boundaryFieldRef()
                            [origPatchi][origPatchFacei] = label(1);
                    }
                }

                break;
            }

            // If this edge is connected to two original faces then add any
            // edge area into non-original connected faces
            case 2:
            {
                forAll(mesh_.edgeFaces()[meshEdgei], edgeFacei)
                {
                    const label facei = mesh_.edgeFaces()[meshEdgei][edgeFacei];

                    const label patchi =
                        mesh_.isInternalFace(facei)
                      ? -1 : pbMesh.patchID()[facei - mesh_.nInternalFaces()];

                    if (patchi != -1 && patchIsOwnerOrig[patchi])
                    {
                        continue;
                    }

                    const label patchFacei =
                        patchi == -1 ? -1 : facei - pbMesh[patchi].start();

                    vector& area =
                        patchi == -1
                      ? SfSf[facei]
                      : SfSf.boundaryFieldRef()[patchi][patchFacei];
                    point& centre =
                        patchi == -1
                      ? CfSf[facei]
                      : CfSf.boundaryFieldRef()[patchi][patchFacei];

                    const label sign =
                        mesh_.faces()[facei].edgeDirection
                        (
                            ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                        );

                    part p(area, centre);
                    p +=
                        sign < 0
                      ? ownerOrigBoundaryEdgeP
                      : -ownerOrigBoundaryEdgeP;

                    area = p.area;
                    centre = p.centre;

                    if (debug && patchi != -1)
                    {
                        tChanged->boundaryFieldRef()
                            [patchi][patchFacei] = label(1);
                    }
                }

                break;
            }

            default:
            {
                FatalErrorInFunction
                    << "Boundary is not manifold"
                    << exit(FatalError);
            }
        }
    }

    // If debugging, check that face changes are in sync
    if (debug)
    {
        const surfaceLabelField::Boundary& changedBf =
            tChanged->boundaryField();

        const surfaceLabelField::Boundary changedNbrBf
        (
            surfaceLabelField::null(),
            changedBf.boundaryNeighbourField()
        );

        bool synchronised = true;

        forAll(changedBf, patchi)
        {
            const fvsPatchLabelField& changedPf = changedBf[patchi];
            const fvsPatchLabelField& changedNbrPf = changedNbrBf[patchi];

            if (!isA<processorFvPatch>(changedPf.patch())) continue;

            forAll(changedPf, patchFacei)
            {
                if (changedPf[patchFacei] != changedNbrPf[patchFacei])
                {
                    const label facei =
                        changedPf.patch().start() + patchFacei;

                    FatalErrorInFunction
                        << "Face not synchronised with centre at "
                        << mesh_.faceCentres()[facei] << endl;

                    synchronised = false;
                }
            }
        }

        if (!synchronised)
        {
            FatalErrorInFunction
                << exit(FatalError);
        }
    }
}


void Foam::fvMeshStitcher::stabiliseOrigPatchFaces
(
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf
) const
{
    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList allOrigPatchIDs = ncb.allOrigPatchIDs();

    forAll(allOrigPatchIDs, i)
    {
        const label origPatchi = allOrigPatchIDs[i];
        const polyPatch& origPp = mesh_.boundaryMesh()[origPatchi];

        forAll(origPp, origPatchFacei)
        {
            part p
            (
                SfBf[origPatchi][origPatchFacei],
                CfBf[origPatchi][origPatchFacei]
            );

            const part smallP
            (
                small*origPp.faceAreas()[origPatchFacei],
                origPp.faceCentres()[origPatchFacei]
            );

            p += smallP;

            SfBf[origPatchi][origPatchFacei] = p.area;
            CfBf[origPatchi][origPatchFacei] = p.centre;
        }
    }
}


void Foam::fvMeshStitcher::intersectNonConformalCyclics
(
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf,
    const bool haveTopology
) const
{
    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList ownerOrigPatchIDs = ncb.ownerOrigPatchIDs();

    // Alias the boundary geometry fields
    surfaceVectorField::Boundary& SfBf = SfSf.boundaryFieldRef();
    surfaceVectorField::Boundary& CfBf = CfSf.boundaryFieldRef();

    // Create storage for and initialise the edge parts of source patches
    List<List<part>> patchEdgeParts(mesh_.boundary().size());
    forAll(ownerOrigPatchIDs, i)
    {
        const label origPatchi = ownerOrigPatchIDs[i];

        patchEdgeParts[origPatchi].resize
        (
            mesh_.boundaryMesh()[origPatchi].nEdges(),
            part(Zero)
        );
    }

    // If topology is specified, also create a boundary field of the original
    // patch indices on the neighbouring interface, as well as the
    // corresponding areas and centres for stabilisation purposes
    tmp<surfaceLabelField::Boundary> tOrigFacesNbrBf;
    tmp<surfaceVectorField::Boundary> tOrigSfNbrBf;
    tmp<surfaceVectorField::Boundary> tOrigCfNbrBf;
    if (haveTopology)
    {
        surfaceLabelField origFaces
        (
            surfaceLabelField::New
            (
                "origFaces",
                mesh_,
                dimensioned<label>(dimless, -1)
            )
        );
        surfaceVectorField origSf("origSf", mesh_.Sf());
        surfaceVectorField origCf("origCf", mesh_.Cf());

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];

            if (!isA<nonConformalFvPatch>(fvp)) continue;

            const nonConformalFvPatch& ncFvp =
                refCast<const nonConformalFvPatch>(fvp);

            origFaces.boundaryFieldRef()[patchi] =
                polyFacesBf[patchi] - ncFvp.origPatch().start();
            origSf.boundaryFieldRef()[patchi] =
                vectorField(mesh_.faceAreas(), polyFacesBf[patchi]);
            origCf.boundaryFieldRef()[patchi] =
                pointField(mesh_.faceCentres(), polyFacesBf[patchi]);
        }

        tOrigFacesNbrBf =
            new surfaceLabelField::Boundary
            (
                surfaceLabelField::null(),
                origFaces.boundaryField().boundaryNeighbourField()
            );
        tOrigSfNbrBf =
            new surfaceVectorField::Boundary
            (
                surfaceVectorField::null(),
                origSf.boundaryField().boundaryNeighbourField()
            );
        tOrigCfNbrBf =
            new surfaceVectorField::Boundary
            (
                surfaceVectorField::null(),
                origCf.boundaryField().boundaryNeighbourField()
            );

        // See note in fvMesh::unconform regarding the transformation
        // of cell centres. The same applies here.
        forAll(tOrigCfNbrBf(), patchi)
        {
            fvsPatchVectorField& Cfp = tOrigCfNbrBf.ref()[patchi];
            if (Cfp.patch().coupled())
            {
                const transformer& t =
                    refCast<const coupledFvPatch>(Cfp.patch())
                   .transform();
                t.invTransform(Cfp, Cfp);
                t.transformPosition(Cfp, Cfp);
            }
        }
    }

    // Do the intersections
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

        const nonConformalCyclicFvPatch& nccFvp =
            refCast<const nonConformalCyclicFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        intersectNonConformalCyclic
        (
            nccFvp,
            polyFacesBf,
            SfBf,
            CfBf,
            tOrigFacesNbrBf,
            tOrigSfNbrBf,
            tOrigCfNbrBf,
            patchEdgeParts[nccFvp.origPatchID()]
        );
    }

    // Construct the boundary edge geometry
    List<part> ownerOrigBoundaryEdgeParts =
        calculateOwnerOrigBoundaryEdgeParts(patchEdgeParts);

    // Add the difference between patch edge parts and all boundary
    // edge parts to the adjacent patch faces. This is an error part.
    forAll(ownerOrigPatchIDs, i)
    {
        const label origPatchi = ownerOrigPatchIDs[i];
        const polyPatch& origPatch = pbMesh[origPatchi];

        const labelList& origPatchEdgeOwnerOrigBoundaryEdges =
            ncb.patchEdgeOwnerOrigBoundaryEdges(origPatchi);

        forAll(patchEdgeParts[origPatchi], origPatchEdgei)
        {
            const label ownerOrigBoundaryEdgei =
                origPatchEdgeOwnerOrigBoundaryEdges[origPatchEdgei];

            part errorP =
                patchEdgeParts[origPatchi][origPatchEdgei];
            errorP -= ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];

            forAll(origPatch.edgeFaces()[origPatchEdgei], patchEdgeFacei)
            {
                const label patchFacei =
                    origPatch.edgeFaces()[origPatchEdgei][patchEdgeFacei];

                part p
                (
                    SfBf[origPatchi][patchFacei],
                    CfBf[origPatchi][patchFacei]
                );
                p += errorP;

                SfBf[origPatchi][patchFacei] = p.area;
                CfBf[origPatchi][patchFacei] = p.centre;
            }
        }
    }

    // Use the boundary edge geometry to correct all edge-connected faces
    applyOwnerOrigBoundaryEdgeParts(SfSf, CfSf, ownerOrigBoundaryEdgeParts);

    // Stabilise
    stabiliseOrigPatchFaces(SfBf, CfBf);
}


template<class NonConformalFvPatch>
inline void Foam::fvMeshStitcher::createNonConformalStabilisationGeometry
(
    const surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf
) const
{
    // Alias the boundary geometry fields
    surfaceVectorField::Boundary& SfBf = SfSf.boundaryFieldRef();
    surfaceVectorField::Boundary& CfBf = CfSf.boundaryFieldRef();

    // Create small stabilisation geometry
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<NonConformalFvPatch>(fvp)) continue;

        const polyPatch& origPp =
            refCast<const nonConformalFvPatch>(fvp).origPatch().patch();

        SfBf[patchi] ==
            vectorField
            (
                small*origPp.faceAreas(),
                polyFacesBf[patchi] - origPp.start()
            );
        CfBf[patchi] ==
            vectorField
            (
                origPp.faceCentres(),
                polyFacesBf[patchi] - origPp.start()
            );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::fvMeshStitcher::geometric() const
{
    bool result = false;

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const fvsPatchScalarField& magSfp =
            mesh_.magSf().boundaryField()[patchi];

        const polyPatch& origPp =
            refCast<const nonConformalFvPatch>(fvp).origPatch().patch();

        const scalarField origMagSfp
        (
            origPp.magFaceAreas(),
            mesh_.polyFacesBf()[patchi] - origPp.start()
        );

        if (max(magSfp/origMagSfp) > rootSmall)
        {
            result = true;
        }
    }

    reduce(result, orOp<bool>());

    return returnReduce(result, orOp<bool>());
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::fvMeshStitcher::openness() const
{
    return
        mag(fvc::surfaceIntegrate(mesh_.Sf()))()*mesh_.V()
       /max
        (
            mag(fvc::surfaceSum(cmptMag(mesh_.Sf())))(),
            (small*sqr(cbrt(mesh_.V())))()
        )();
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::fvMeshStitcher::volumeConservationError(const label n) const
{
    if (0 > n || n > 1)
    {
        FatalErrorInFunction
            << "Can only compute volume conservation error for this time, or "
            << "the previous time" << exit(FatalError);
    }

    const surfaceScalarField& phi = mesh_.phi().oldTime(n);

    const dimensionedScalar deltaT =
        n == 0 ? mesh_.time().deltaT() : mesh_.time().deltaT0();

    const volScalarField::Internal& V = n == 0 ? mesh_.V() : mesh_.V0();
    const volScalarField::Internal& V0 = n == 0 ? mesh_.V0() : mesh_.V00();

    return fvc::surfaceIntegrate(phi*deltaT)() - (V - V0)/mesh_.V();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshStitcher::fvMeshStitcher(fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshStitcher::~fvMeshStitcher()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshStitcher::stitches() const
{
    forAll(mesh_.boundary(), patchi)
    {
        if (isA<nonConformalFvPatch>(mesh_.boundary()[patchi]))
        {
            return true;
        }
    }

    return false;
}


void Foam::fvMeshStitcher::updateMesh(const polyTopoChangeMap&)
{}


void Foam::fvMeshStitcher::movePoints()
{}


bool Foam::fvMeshStitcher::disconnect
(
    const bool changing,
    const bool geometric
)
{
    if (!stitches() || (changing && !mesh_.dynamic()))
    {
        return false;
    }

    // Is this mesh coupled?
    const bool coupled = Pstream::parRun() || !mesh_.time().processorCase();

    if (coupled && geometric)
    {
        Info<< indent << typeName << ": Disconnecting" << incrIndent << endl;
    }

    if (changing)
    {
        // Pre-conform surface fields. This splits the original and cyclic
        // parts of the interface fields into separate boundary fields, with
        // both sets of values store on the original faces. The original field
        // overwrites the existing boundary values, whilst the cyclic field is
        // stored as a separate field for use later.
        preConformSurfaceFields();
    }

    // Conform the mesh
    if (mesh_.moving())
    {
        surfaceScalarField phi(mesh_.phi());
        conformCorrectMeshPhi(phi);
        mesh_.conform(phi);
    }
    else
    {
        mesh_.conform();
    }

    // Resize all the affected patch fields
    resizePatchFields<VolField>();
    resizePatchFields<SurfaceField>();

    // Prevent hangs caused by processor cyclic patches using mesh geometry
    mesh_.deltaCoeffs();

    if (coupled && geometric)
    {
        const volScalarField::Internal o(openness());
        Info<< indent << "Cell min/avg/max openness = "
            << gMin(o) << '/' << gAverage(o) << '/' << gMax(o) << endl;

        for (label i = 0; mesh_.moving() && i <= mesh_.phi().nOldTimes(); ++ i)
        {
            const volScalarField::Internal vce(volumeConservationError(i));
            Info<< indent << "Cell min/avg/max ";
            for (label j = 0; j < i; ++ j) Info<< "old-";
            Info<< (i ? "time " : "") << "volume conservation error = "
                << gMin(vce) << '/' << gAverage(vce) << '/' << gMax(vce)
                << endl;
        }

    }

    if (coupled && geometric)
    {
        Info<< decrIndent;
    }

    // Create null polyTopoChangeMap
    const polyTopoChangeMap map(mesh_);

    meshObjects::topoChange<fvMesh>(mesh_, map);
    meshObjects::topoChange<lduMesh>(mesh_, map);

    const_cast<Time&>(mesh_.time()).functionObjects().topoChange(map);

    return true;
}


bool Foam::fvMeshStitcher::connect
(
    const bool changing,
    const bool geometric,
    const bool load
)
{
    if (!stitches() || (changing && !mesh_.dynamic()))
    {
        return false;
    }

    // Is this mesh coupled?
    const bool coupled = Pstream::parRun() || !mesh_.time().processorCase();

    // Create a copy of the conformal poly face addressing
    surfaceLabelField::Boundary polyFacesBf
    (
        surfaceLabelField::null(),
        mesh_.polyFacesBf()
    );

    // If starting up then load topology from disk, if it is available
    bool haveTopology = false;
    if (load)
    {
        const IOobject polyFacesBfIO(mesh_.polyFacesBfIO(IOobject::MUST_READ));

        if (fileHandler().isFile(polyFacesBfIO.objectPath(false)))
        {
            haveTopology = true;

            // Read the boundary field but then set default values for conformal
            // patches as these patches will have had a uniform invalid index
            // set in order to save disk space
            surfaceLabelField::Boundary polyFacesBfRead
            (
                surfaceLabelField::null(),
                surfaceLabelField(polyFacesBfIO, mesh_).boundaryField()
            );
            forAll(mesh_.boundary(), patchi)
            {
                if (isA<nonConformalFvPatch>(mesh_.boundary()[patchi]))
                {
                    polyFacesBf[patchi] = polyFacesBfRead[patchi];
                }
            }
        }
    }

    // Access all the intersections in advance. Makes the log nicer.
    if (coupled && (geometric || !haveTopology))
    {
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];

            if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

            const nonConformalCyclicFvPatch& nccFvp =
                refCast<const nonConformalCyclicFvPatch>(fvp);

            if (!nccFvp.owner()) continue;

            nccFvp.nonConformalCyclicPatch().intersection();
        }
    }

    if (coupled && (geometric || !haveTopology))
    {
        Info<< indent << typeName << ": Connecting" << incrIndent << endl;
    }

    // Create copies of geometry fields to be modified
    surfaceVectorField Sf(mesh_.Sf().cloneUnSliced()());
    surfaceVectorField Cf(mesh_.Cf().cloneUnSliced()());

    // Construct non-conformal geometry
    if (coupled && (geometric || !haveTopology))
    {
        // Do the intersection and create the non-conformal cyclic faces
        intersectNonConformalCyclics(polyFacesBf, Sf, Cf, haveTopology);

        // Create stabilisation geometry on any error faces specified in the
        // addressing
        createNonConformalStabilisationGeometry<nonConformalErrorFvPatch>
        (
            polyFacesBf,
            Sf,
            Cf
        );
    }
    else
    {
        // If not coupled, then create stabilisation geometry on all
        // non-conformal faces specified in the addressing
        createNonConformalStabilisationGeometry<nonConformalFvPatch>
        (
            polyFacesBf,
            Sf,
            Cf
        );
    }

    // If the mesh is moving then create any additional non-conformal geometry
    // necessary to correct the mesh fluxes
    if (mesh_.moving())
    {
        createNonConformalCorrectMeshPhiGeometry(polyFacesBf, Sf, Cf);
    }

    // Make the mesh non-conformal
    if (mesh_.moving())
    {
        surfaceScalarField phi(mesh_.phi());
        unconformCorrectMeshPhi(polyFacesBf, Sf, Cf, phi);
        mesh_.unconform(polyFacesBf, Sf, Cf, phi);
    }
    else
    {
        mesh_.unconform(polyFacesBf, Sf, Cf);
    }

    // Resize all the affected patch fields
    resizePatchFields<VolField>();
    resizePatchFields<SurfaceField>();

    if (changing)
    {
        // Post-non-conform surface fields. This reconstructs the original and
        // cyclic parts of the interface fields from separate original and
        // cyclic parts. The original part was store in the same field, whilst
        // the cyclic part was separately registered.
        postNonConformSurfaceFields();

        // Volume fields are assumed to be intensive. So, the value on a face
        // which has changed in size can be retained without modification. New
        // faces need values to be set. This is done by evaluating all the
        // nonConformalCoupled patch fields.
        evaluateVolFields();

        // Do special post-non-conformation for surface velocities.
        postNonConformSurfaceVelocities();
    }

    // Prevent hangs caused by processor cyclic patches using mesh geometry
    mesh_.deltaCoeffs();

    if (coupled && geometric)
    {
        const volScalarField::Internal o(openness());
        Info<< indent << "Cell min/avg/max openness = "
            << gMin(o) << '/' << gAverage(o) << '/' << gMax(o) << endl;

        for (label i = 0; mesh_.moving() && i <= mesh_.phi().nOldTimes(); ++ i)
        {
            const volScalarField::Internal vce(volumeConservationError(i));
            Info<< indent << "Cell min/avg/max ";
            for (label j = 0; j < i; ++ j) Info<< "old-";
            Info<< (i ? "time " : "") << "volume conservation error = "
                << gMin(vce) << '/' << gAverage(vce) << '/' << gMax(vce)
                << endl;
        }

        if (mesh_.moving())
        {
            // Create a boundary field of the imbalance between the mesh fluxes
            // on either side of interfaces (like phiErrorb above)
            surfaceScalarField::Boundary mfe
            (
                surfaceScalarField::Internal::null(),
                mesh_.phi().boundaryField()
            );
            mfe += mesh_.phi().boundaryField().boundaryNeighbourField();

            // Determine the number of non-processor patches
            label nNonProcPatches = 0;
            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                if (!isA<processorFvPatch>(fvp))
                {
                    if (nNonProcPatches != patchi)
                    {
                        FatalErrorInFunction
                            << "Processor patches do not follow non-processor "
                            << "patches" << exit(FatalError);
                    }

                    nNonProcPatches ++;
                }
            }

            // Work out the min, avg and max error for all non-conformal
            // cyclics, including any errors on referred processor cyclics
            scalarList minMfe(nNonProcPatches, vGreat);
            scalarField sumMfe(nNonProcPatches, 0), nSumMfe(nNonProcPatches, 0);
            scalarList maxMfe(nNonProcPatches, -vGreat);
            forAll(mfe, patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                const label nccPatchi =
                    isA<nonConformalCyclicFvPatch>(fvp)
                  ? refCast<const nonConformalCyclicFvPatch>(fvp)
                   .index()
                  : isA<nonConformalProcessorCyclicFvPatch>(fvp)
                  ? refCast<const nonConformalProcessorCyclicFvPatch>(fvp)
                   .referPatchID()
                  : -1;

                if (nccPatchi != -1)
                {
                    minMfe[nccPatchi] =
                        min(minMfe[nccPatchi], min(mfe[patchi]));
                    sumMfe[nccPatchi] += sum(mfe[patchi]);
                    nSumMfe[nccPatchi] += mfe[patchi].size();
                    maxMfe[nccPatchi] =
                        max(maxMfe[nccPatchi], max(mfe[patchi]));
                }
            }
            reduce(minMfe, ListOp<minOp<scalar>>());
            reduce(sumMfe, ListOp<minOp<scalar>>());
            reduce(nSumMfe, ListOp<minOp<scalar>>());
            reduce(maxMfe, ListOp<maxOp<scalar>>());

            // Report
            forAll(minMfe, patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                if
                (
                    isA<nonConformalCyclicFvPatch>(fvp)
                 && refCast<const nonConformalCyclicFvPatch>(fvp).owner()
                 && nSumMfe[patchi] > 0
                )
                {
                    Info<< indent << fvp.name()
                        << " min/avg/max mesh flux error = " << minMfe[patchi]
                        << '/' << sumMfe[patchi]/(nSumMfe[patchi] + vSmall)
                        << '/' << maxMfe[patchi] << endl;
                }
            }
        }
    }

    if (coupled && (geometric || !haveTopology))
    {
        Info<< decrIndent;
    }

    // Create null polyTopoChangeMap
    const polyTopoChangeMap map(mesh_);

    meshObjects::topoChange<fvMesh>(mesh_, map);
    meshObjects::topoChange<lduMesh>(mesh_, map);

    const_cast<Time&>(mesh_.time()).functionObjects().topoChange(map);

    return true;
}


void Foam::fvMeshStitcher::reconnect(const bool geometric) const
{
    if (mesh_.conformal() || geometric == this->geometric())
    {
        return;
    }

    // Create a copy of the conformal poly face addressing
    surfaceLabelField::Boundary polyFacesBf
    (
        surfaceLabelField::null(),
        mesh_.polyFacesBf()
    );

    // Undo all non-conformal changes and clear all geometry and topology
    mesh_.conform();

    // Create copies of geometry fields to be modified
    surfaceVectorField Sf(mesh_.Sf().cloneUnSliced()());
    surfaceVectorField Cf(mesh_.Cf().cloneUnSliced()());

    // Construct non-conformal geometry
    if (geometric)
    {
        // Do the intersection and create the non-conformal cyclic faces
        intersectNonConformalCyclics(polyFacesBf, Sf, Cf, true);

        // Create stabilisation geometry on any error faces specified in the
        // addressing
        createNonConformalStabilisationGeometry<nonConformalErrorFvPatch>
        (
            polyFacesBf,
            Sf,
            Cf
        );
    }
    else
    {
        // If not coupled, then create stabilisation geometry on all
        // non-conformal faces specified in the addressing
        createNonConformalStabilisationGeometry<nonConformalFvPatch>
        (
            polyFacesBf,
            Sf,
            Cf
        );
    }

    // Apply changes to the mesh
    mesh_.unconform(polyFacesBf, Sf, Cf);

    // Prevent hangs caused by processor cyclic patches using mesh geometry
    mesh_.deltaCoeffs();

    // Create null polyTopoChangeMap
    const polyTopoChangeMap map(mesh_);

    meshObjects::topoChange<fvMesh>(mesh_, map);
    meshObjects::topoChange<lduMesh>(mesh_, map);

    const_cast<Time&>(mesh_.time()).functionObjects().topoChange(map);
}


// ************************************************************************* //
