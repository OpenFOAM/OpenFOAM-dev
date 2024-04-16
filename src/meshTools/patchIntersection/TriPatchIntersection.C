/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "TriPatchIntersection.H"
#include "barycentricTensor2D.H"
#include "boundSphere.H"
#include "cpuTime.H"
#include "indexedOctree.H"
#include "OFstream.H"
#include "treeDataPrimitivePatch.H"
#include "triIntersect.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::checkPatchFace
(
    const label patchFacei,
    const bool isSrc
) const
{
    #ifdef FULLDEBUG

    const labelList& patchFaceTris =
        isSrc ? srcFaceTris_[patchFacei] : tgtFaceTris_[patchFacei];

    forAll(patchFaceTris, patchFaceTrii)
    {
        const label trii = patchFaceTris[patchFaceTrii];

        if (triPoints_[trii] != FixedList<label, 3>({-1, -1, -1}))
        {
            forAll(triEdges_[trii], triEdgei)
            {
                const label edgei = triEdges_[trii][triEdgei];

                // Check that there are no duplicate points
                const edge e = triEdgePoints(trii, triEdgei);
                if (e[0] == e[1])
                {
                    FatalErrorInFunction
                        << "Tri #" << trii << "'s tri-edge #" << triEdgei
                        << " (edge #" << triEdges_[trii][triEdgei] << ") has "
                        << "duplicate points " << e << exit(FatalError);
                }

                // Check that this edge also references this tri
                const label edgeTrii = findIndex(edgeTris_[edgei], trii);
                if (edgeTrii == -1)
                {
                    FatalErrorInFunction
                        << "Tri #" << trii << " references edge #" << edgei
                        << " but the reverse is not true" << exit(FatalError);
                }

                // If there is a connected tri ...
                const label otherTrii = edgeTris_[edgei][!edgeTrii];
                if (otherTrii != -1)
                {
                    // Check that the connected tri also references this edge
                    const label otherTriEdgei =
                        findIndex(triEdges_[otherTrii], edgei);
                    if (otherTriEdgei == -1)
                    {
                        FatalErrorInFunction
                            << "Edge #" << edgei << " references tri #"
                            << otherTrii << " but the reverse is not true"
                            << exit(FatalError);
                    }

                    // Check that the edge is correspondingly oriented in its
                    // two connected tris
                    const edge otherE = triEdgePoints(otherTrii, otherTriEdgei);
                    if (edge::compare(e, otherE) != -1)
                    {
                        FatalErrorInFunction
                            << "Edge #" << edgei << ' ' << edgePoints(edgei)
                            << " is not the same in adjacent tri #" << trii
                            << ' ' << triPoints(trii) << " and tri #"
                            << otherTrii << ' ' << triPoints(otherTrii)
                            << exit(FatalError);
                    }
                }

                // Check edge-points and point-edges ... !!!

                // If this is a front edge ...
                if (edgeFrontEdges_[edgei] != -1)
                {
                    // Check that the front arrays are consistent
                    const label frontEdgei = edgeFrontEdges_[edgei];
                    if
                    (
                        frontEdgeEdges_.size() <= frontEdgei
                     || frontEdgeEdges_[frontEdgei] != edgei
                    )
                    {
                        FatalErrorInFunction
                            << "Edge #" << edgei << " is marked as part of the "
                            << "front but is not in the front edge list"
                            << exit(FatalError);
                    }

                    // Check that the edge connects source to target
                    const label trii0 = edgeTris_[edgei][0];
                    const label trii1 = edgeTris_[edgei][1];
                    if
                    (
                        trii0 == -1
                     || trii1 == -1
                     || (triSrcFace_[trii0] == -1) == (triSrcFace_[trii1] == -1)
                    )
                    {
                        FatalErrorInFunction
                            << "Front edge #" << edgei
                            << " does not connect a source-tri to a target-tri"
                            << exit(FatalError);
                    }
                }

                // If this is a star edge ... !!!
            }

            // If this is a candidate tri ... !!!

            // If this is a marked tri ... !!!
        }
    }

    #endif
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::checkPatchEdge
(
    const label patchEdgei,
    const bool isSrc
) const
{
    #ifdef FULLDEBUG

    const labelListList& patchEdgePatchFaces =
        isSrc ? this->srcPatch_.edgeFaces() : this->tgtPatch_.edgeFaces();

    const labelList patchFaceis = patchEdgePatchFaces[patchEdgei];

    forAll(patchFaceis, i)
    {
        checkPatchFace(patchFaceis[i], isSrc);
    }

    #endif
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::checkPatchFaces
(
    const bool isSrc
) const
{
    #ifdef FULLDEBUG

    const List<DynamicList<label>>& patchFaceTris =
        isSrc ? srcFaceTris_ : tgtFaceTris_;

    forAll(patchFaceTris, patchFacei)
    {
        checkPatchFace(patchFacei, isSrc);
    }

    #endif
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::removeEdge
(
    const label edgei
)
{
    #ifdef FULLDEBUG
    if
    (
        edgeTris_[edgei] != labelPair(-1, -1)
     && intersectEdgeFaces_[edgei] != labelPair(-1, -1)
    )
    {
        FatalErrorInFunction
            << "Attempted to remove edge #" << edgei << " which is still "
            << "connected to triangles " << edgeTris_[edgei] << " and faces "
            << intersectEdgeFaces_[edgei]
            << exit(FatalError);
    }
    #endif

    removedEdges_.append(edgei);

    if (edgeFrontEdges_[edgei] != -1)
    {
        frontEdgeEdges_[edgeFrontEdges_[edgei]] = -1;
        edgeFrontEdges_[edgei] = -1;
    }
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::removeTri
(
    const label trii
)
{
    triPoints_[trii] = triFace(-1, -1, -1);

    {
        const bool isSrc = triSrcFace_[trii] != -1;

        const label patchFacei = isSrc ? triSrcFace_[trii] : triTgtFace_[trii];

        DynamicList<label>& patchFaceTris =
            isSrc ? srcFaceTris_[patchFacei] : tgtFaceTris_[patchFacei];

        label patchFaceTrii = patchFaceTris.size() - 1;

        for (; patchFaceTris[patchFaceTrii] != trii; -- patchFaceTrii);

        for (; patchFaceTrii < patchFaceTris.size() - 1; ++ patchFaceTrii)
        {
            patchFaceTris[patchFaceTrii] = patchFaceTris[patchFaceTrii + 1];
        }

        patchFaceTris.resize(patchFaceTris.size() - 1);
    }

    triSrcFace_[trii] = -1;
    triTgtFace_[trii] = -1;

    if (triCandidateTris_[trii] != -1)
    {
        candidateTriTris_[triCandidateTris_[trii]] = -1;
        triCandidateTris_[trii] = -1;
    }

    if (triMarkedTris_[trii] != -1)
    {
        markedTriTris_[triMarkedTris_[trii]] = -1;
        triMarkedTris_[trii] = -1;
    }

    forAll(triEdges_[trii], i)
    {
        const label edgei = triEdges_[trii][i];

        edgeTris_[edgei][edgeTris_[edgei][1] == trii] = -1;

        if
        (
            edgeTris_[edgei] == labelPair(-1, -1)
         && intersectEdgeFaces_[edgei] == labelPair(-1, -1)
        )
        {
            removeEdge(edgei);
        }
    }

    triEdges_[trii] = FixedList<label, 3>({-1, -1, -1});

    removedTris_.append(trii);
}


template<class SrcPatchType, class TgtPatchType>
Foam::label Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::newTrii()
{
    if (removedTris_.size())
    {
        return removedTris_.remove();
    }
    else
    {
        triPoints_.append(triFace(-1, -1, -1));
        triEdges_.append(FixedList<label, 3>({-1, -1, -1}));
        triSrcFace_.append(-1);
        triTgtFace_.append(-1);
        triCandidateTris_.append(-1);
        triMarkedTris_.append(-1);
        return triPoints_.size() - 1;
    }
}


template<class SrcPatchType, class TgtPatchType>
Foam::label Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::newEdgei()
{
    if (removedEdges_.size())
    {
        return removedEdges_.remove();
    }
    else
    {
        edgeTris_.append({-1, -1});
        intersectEdgeFaces_.append({-1, -1});
        edgeFrontEdges_.append(-1);
        return edgeTris_.size() - 1;
    }
}


template<class SrcPatchType, class TgtPatchType>
Foam::label Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::newPointi()
{
    const label pointi = srcPoints_.size();

    srcPoints_.append(point::uniform(NaN));
    srcPointNormals_.append(vector::uniform(NaN));
    tgtPoints_.append(point::uniform(NaN));
    pointPoints_.append(pointi);
    this->pointSrcFaces_.append(-1);
    this->pointTgtFaces_.append(-1);
    this->pointSrcEdges_.append(-1);
    this->pointTgtEdges_.append(-1);
    this->pointSrcPoints_.append(-1);
    this->pointTgtPoints_.append(-1);
    return pointi;
}


template<class SrcPatchType, class TgtPatchType>
Foam::label Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triPoint
(
    const label trii,
    const label triPointi
) const
{
    return pointPoints_[triPoints_[trii][triPointi]];
}


template<class SrcPatchType, class TgtPatchType>
Foam::triFace Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triPoints
(
    const label trii
) const
{
    triFace result;
    forAll(triPoints_[trii], triPointi)
    {
        result[triPointi] = triPoint(trii, triPointi);
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
template<class Type>
Foam::FixedList<Type, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triPointValues
(
    const label trii,
    const UList<Type> values
) const
{
    FixedList<Type, 3> result;
    forAll(result, triPointi)
    {
        result[triPointi] = values[triPoint(trii, triPointi)];
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::FixedList<bool, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triOwns
(
    const label trii
) const
{
    FixedList<bool, 3> result;
    forAll(triEdges_[trii], triEdgei)
    {
        result[triEdgei] = edgeTris_[triEdges_[trii][triEdgei]][0] == trii;
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::FixedList<Foam::label, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triOtherTriPoints
(
    const label trii,
    const label otherTrii
) const
{
    FixedList<label, 3> result({-1, -1, -1});
    forAll(triPoints_[trii], triPointi)
    {
        forAll(triPoints_[otherTrii], otherTriPointi)
        {
            if
            (
                triPoint(trii, triPointi)
             == triPoint(otherTrii, otherTriPointi)
            )
            {
                result[triPointi] = otherTriPointi;
            }
        }
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::FixedList<Foam::label, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triOtherTriEdges
(
    const label trii,
    const label otherTrii
) const
{
    FixedList<label, 3> result({-1, -1, -1});
    forAll(triEdges_[trii], triEdgei)
    {
        forAll(triEdges_[otherTrii], otherTriEdgei)
        {
            if
            (
                triEdges_[trii][triEdgei]
             == triEdges_[otherTrii][otherTriEdgei]
            )
            {
                result[triEdgei] = otherTriEdgei;
            }
        }
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::edge Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::triEdgePoints
(
    const label trii,
    const label triEdgei
) const
{
    return edge
    (
        triPoint(trii, triEdgei),
        triPoint(trii, (triEdgei + 1) % 3)
    );
}


template<class SrcPatchType, class TgtPatchType>
Foam::edge Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::edgePoints
(
    const label edgei
) const
{
    const label edgeSidei = edgeTris_[edgei][0] == -1;

    const label trii = edgeTris_[edgei][edgeSidei];
    const label triEdgei = findIndex(triEdges_[trii], edgei);

    const edge e = triEdgePoints(trii, triEdgei);

    return edgeSidei == 0 ? e : e.reverseEdge();
}


template<class SrcPatchType, class TgtPatchType>
Foam::label
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePoint
(
    const label patchFacei,
    const label patchFacePointi,
    const bool isSrc
) const
{
    return
        isSrc
      ? this->srcPointPoints_
        [
            this->srcPatch_.localFaces()[patchFacei][patchFacePointi]
        ]
      : this->tgtPointPoints_
        [
            this->tgtPatch_.localFaces()[patchFacei][patchFacePointi]
        ];
}


template<class SrcPatchType, class TgtPatchType>
Foam::triFace
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePoints
(
    const label patchFacei,
    const bool isSrc
) const
{
    triFace result;
    forAll(result, patchFacePointi)
    {
        result[patchFacePointi] =
            patchFacePoint(patchFacei, patchFacePointi, isSrc);
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
template<class Type>
Foam::FixedList<Type, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePointValues
(
    const label patchFacei,
    const bool isSrc,
    const UList<Type>& values
) const
{
    FixedList<vector, 3> result;
    forAll(result, patchFacePointi)
    {
        result[patchFacePointi] =
            values[patchFacePoint(patchFacei, patchFacePointi, isSrc)];
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::FixedList<bool, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFaceOwns
(
    const label patchFacei,
    const bool isSrc
) const
{
    const UList<triFace>& localFaces =
        isSrc ? this->srcPatch_.localFaces() : this->tgtPatch_.localFaces();
    const labelListList& faceEdges =
        isSrc ? this->srcPatch_.faceEdges() : this->tgtPatch_.faceEdges();
    const edgeList& localEdges =
        isSrc ? this->srcPatch_.edges() : this->tgtPatch_.edges();
    FixedList<bool, 3> result;
    forAll(result, i)
    {
        const edge& e = localEdges[faceEdges[patchFacei][i]];
        const edge fe = localFaces[patchFacei].faceEdge(i);
        result[i] = edge::compare(e, fe) > 0;
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::FixedList<Foam::label, 3>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::
patchFaceOtherPatchFacePoints
(
    const label patchFacei,
    const label otherPatchFacei,
    const bool isSrc
) const
{
    const triFace pointis = patchFacePoints(patchFacei, isSrc);
    const triFace otherPointis = patchFacePoints(otherPatchFacei, !isSrc);
    FixedList<label, 3> result({-1, -1, -1});
    forAll(pointis, i)
    {
        forAll(otherPointis, otheri)
        {
            if (pointis[i] == otherPointis[otheri])
            {
                result[i] = otheri;
            }
        }
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::label
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePatchPoint
(
    const label patchFacei,
    const label patchFacePointi,
    const bool isSrc
) const
{
    return
        isSrc
      ? this->srcPatch_.localFaces()[patchFacei][patchFacePointi]
      : this->tgtPatch_.localFaces()[patchFacei][patchFacePointi];
}


template<class SrcPatchType, class TgtPatchType>
Foam::triFace
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePatchPoints
(
    const label patchFacei,
    const bool isSrc
) const
{
    triFace result;
    forAll(result, i)
    {
        result[i] = patchFacePatchPoint(patchFacei, i, isSrc);
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::label
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePatchEdge
(
    const label patchFacei,
    const label patchFaceEdgei,
    const bool isSrc
) const
{
    return
        isSrc
      ? this->srcPatch_.faceEdges()[patchFacei][patchFaceEdgei]
      : this->tgtPatch_.faceEdges()[patchFacei][patchFaceEdgei];
}


template<class SrcPatchType, class TgtPatchType>
Foam::triFace
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::patchFacePatchEdges
(
    const label patchFacei,
    const bool isSrc
) const
{
    triFace result;
    forAll(result, i)
    {
        result[i] = patchFacePatchEdge(patchFacei, i, isSrc);
    }
    return result;
}


template<class SrcPatchType, class TgtPatchType>
Foam::label
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::edgePatchEdge
(
    const label edgei,
    const bool isSrc
) const
{
    const labelListList& patchPointPatchEdges =
        isSrc ? this->srcPatch_.pointEdges() : this->tgtPatch_.pointEdges();

    const labelList& pointPatchPoints =
        isSrc ? this->pointSrcPoints_ : this->pointTgtPoints_;
    const labelList& pointPatchEdges =
        isSrc ? this->pointSrcEdges_ : this->pointTgtEdges_;

    auto nPointPatchEdges = [&](const label pointi)
    {
        const label patchPointi = pointPatchPoints[pointi];

        if (patchPointi == -1)
        {
            return pointPatchEdges[pointi] == -1 ? label(0) : label(1);
        }
        else
        {
            return patchPointPatchEdges[patchPointi].size();
        }
    };

    auto pointPatchEdge = [&]
    (
        const label pointi,
        const label i
    )
    {
        const label patchPointi = pointPatchPoints[pointi];

        if (patchPointi == -1)
        {
            return pointPatchEdges[pointi];
        }
        else
        {
            return patchPointPatchEdges[patchPointi][i];
        }
    };

    const edge& e = edgePoints(edgei);
    const label pointi0 = e.first(), pointi1 = e.last();

    label patchEdgei = -1;

    for (label i0 = 0; i0 < nPointPatchEdges(pointi0); ++ i0)
    {
        for (label i1 = 0; i1 < nPointPatchEdges(pointi1); ++ i1)
        {
            const label patchEdgei0 = pointPatchEdge(pointi0, i0);
            const label patchEdgei1 = pointPatchEdge(pointi1, i1);

            if (patchEdgei0 == patchEdgei1)
            {
                if (patchEdgei != -1 && patchEdgei != patchEdgei0)
                {
                    FatalErrorInFunction
                        << "Edge #" << edgei << " (points " << e
                        << ") is associated with two or more different "
                        << (isSrc ? "source" : "target")
                        << " patch edges. This should not be possible."
                        << exit(FatalError);
                }

                patchEdgei = patchEdgei0;
            }
        }
    }

    return patchEdgei;
}


template<class SrcPatchType, class TgtPatchType>
Foam::labelPair
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::edgePatchEdges
(
    const label edgei
) const
{
    return labelPair(edgePatchEdge(edgei, true), edgePatchEdge(edgei, false));
}


template<class SrcPatchType, class TgtPatchType>
Foam::label Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::addTri
(
    const triFace& pointis,
    const FixedList<label, 3>& edgeis,
    const label patchFacei,
    const bool isSrc
)
{
    // Check that edge points match the existing mesh
    #ifdef FULLDEBUG
    forAll(edgeis, triEdgei)
    {
        const label edgei = edgeis[triEdgei];

        if (edgei == -1) continue;

        const label otherTrii = edgeTris_[edgei][edgeTris_[edgei][0] == -1];

        if (otherTrii == -1) continue;

        const label otherTriEdgei = findIndex(triEdges_[otherTrii], edgei);

        const edge eThis = pointis.faceEdge(triEdgei);
        const edge eOther = triEdgePoints(otherTrii, otherTriEdgei);

        if (edge::compare(eThis, eOther) != -1)
        {
            FatalErrorInFunction
                << "Edge #" << edgei << ' ' << edgePoints(edgei)
                << " is not the same in the new tri " << pointis << ' '
                << edgeis << " and the existing adjacent tri "
                << triPoints(otherTrii) << ' ' << triEdges_[otherTrii]
                << exit(FatalError);
        }
    }
    #endif

    // Get the new triangle label
    const label trii = newTrii();

    // Set the points
    triPoints_[trii] = pointis;

    // Set the edges, creating new ones where specified, and completing the
    // edge-tri associations
    triEdges_[trii] = edgeis;
    forAll(triEdges_[trii], triEdgei)
    {
        label& edgei = triEdges_[trii][triEdgei];

        if (edgei == -1)
        {
            edgei = newEdgei();
        }

        edgeTris_[edgei][edgeTris_[edgei][0] != -1] = trii;
    }

    // Set the patch face association
    if (isSrc)
    {
        triSrcFace_[trii] = patchFacei;
        srcFaceTris_[patchFacei].append(trii);
    }
    else
    {
        triTgtFace_[trii] = patchFacei;
        tgtFaceTris_[patchFacei].append(trii);
    }

    return trii;
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::flipEdge
(
    const label edgei
)
{
    // Return if this edge does not have two adjacent triangles
    forAll(edgeTris_[edgei], edgeTrii)
    {
        const label trii = edgeTris_[edgei][edgeTrii];

        if (trii == -1) return;
    }

    // Check that the flip doesn't break patch face associations
    #ifdef FULLDEBUG
    if (edgePatchEdges(edgei) != labelPair(-1, -1))
    {
        FatalErrorInFunction
            << "Flipping an original edge"
            << exit(FatalError);
    }
    #endif

    // Store stuff that will disappear when the existing triangles are modified
    label pointi0 = -1;
    label pointi1 = -1;
    labelPair pointiOpps(-1, -1);
    labelPair edgei0s(-1, -1), trii0s(-1, -1);
    labelPair edgei1s(-1, -1), trii1s(-1, -1);
    forAll(edgeTris_[edgei], edgeTrii)
    {
        const label trii = edgeTris_[edgei][edgeTrii];
        const label triEdgei = findIndex(triEdges_[trii], edgei);

        pointi0 = triPoint(trii, (triEdgei + edgeTrii) % 3);
        pointi1 = triPoint(trii, (triEdgei + !edgeTrii) % 3);
        pointiOpps[edgeTrii] = triPoint(trii, (triEdgei + 2) % 3);

        edgei0s[edgeTrii] = triEdges_[trii][(triEdgei + !edgeTrii + 1) % 3];
        trii0s[edgeTrii] =
            edgeTris_
            [edgei0s[edgeTrii]]
            [edgeTris_[edgei0s[edgeTrii]][0] == trii];
        edgei1s[edgeTrii] = triEdges_[trii][(triEdgei + edgeTrii + 1) % 3];
        trii1s[edgeTrii] =
            edgeTris_
            [edgei1s[edgeTrii]]
            [edgeTris_[edgei1s[edgeTrii]][0] == trii];
    }

    // Modify the triangles
    {
        const label trii0 = edgeTris_[edgei][0], trii1 = edgeTris_[edgei][1];

        triPoints_[trii0] = {pointi0, pointiOpps[1], pointiOpps[0]};
        triEdges_[trii0] = {edgei0s[1], edgei, edgei0s[0]};
        edgeTris_[edgei0s[1]][edgeTris_[edgei0s[1]][1] == trii1] = trii0;

        triPoints_[trii1] = {pointi1, pointiOpps[0], pointiOpps[1]};
        triEdges_[trii1] = {edgei1s[0], edgei, edgei1s[1]};
        edgeTris_[edgei1s[0]][edgeTris_[edgei1s[0]][1] == trii0] = trii1;
    }
}


template<class SrcPatchType, class TgtPatchType>
Foam::scalar
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::circumDistSqr
(
    const label trii,
    const label pointi
) const
{
    // !!! It may be preferable to do this in the plane of the patch face,
    // rather than the current triangle, but that would mean implementing
    // custom circumcircle methods.

    const pointField& points =
        triSrcFace_[trii] != -1 ? srcPoints_ : tgtPoints_;

    const triPointRef t = triPoints(trii).tri(points);

    const Tuple2<point, scalar> circle = t.circumCircle();

    return magSqr(points[pointi] - circle.first()) - sqr(circle.second());
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::insertPoints
(
    const label insertionTriOrEdgei,
    const bool isTri,
    const UList<label>& pointis,
    UList<label>& insertionEdgeis,
    const UList<label>& fixedEdgeis
)
{
    // Clear the insertion edge list
    insertionEdgeis = -1;

    // If inserting into a tri then initialise the candidates with the given
    // insertion tri
    if (isTri)
    {
        candidateTriTris_.append(insertionTriOrEdgei);
        triCandidateTris_[insertionTriOrEdgei] = 0;
    }

    // If inserting into an edge, get the insertion edge label
    if (!isTri)
    {
        insertionEdgeis[0] = insertionTriOrEdgei;
    }

    // Check that we are inserting unambiguously into the source or target side
    #ifdef FULLDEBUG
    if (!isTri)
    {
        const label trii0 = edgeTris_[insertionTriOrEdgei][0];
        const label trii1 = edgeTris_[insertionTriOrEdgei][1];

        if
        (
            trii0 != -1
         && trii1 != -1
         && (triSrcFace_[trii0] == -1) != (triSrcFace_[trii1] == -1)
        )
        {
            FatalErrorInFunction
                << "Attempted insertion into front edge #"
                << insertionTriOrEdgei << exit(FatalError);
        }
    }
    #endif

    // Determine the patch associations for the insertion tri or edge
    const label adjacentTri =
        isTri
      ? insertionTriOrEdgei
      : edgeTris_[insertionTriOrEdgei][edgeTris_[insertionTriOrEdgei][0] == -1];
    const bool isSrc = triSrcFace_[adjacentTri] != -1;
    const label patchEdgei =
        isTri
      ? -1
      : edgePatchEdge(insertionTriOrEdgei, isSrc);
    const label patchFacei =
        patchEdgei != -1
      ? -1
      : (isSrc ? triSrcFace_[adjacentTri] : triTgtFace_[adjacentTri]);

    // Get connectivity arrays that will be modified during insertion
    List<DynamicList<label>>& patchEdgePoints =
        isSrc ? this->srcEdgePoints_ : this->tgtEdgePoints_;
    DynamicList<label>& pointPatchEdges =
        isSrc ? this->pointSrcEdges_ : this->pointTgtEdges_;
    DynamicList<label>& pointPatchFaces =
        isSrc ? this->pointSrcFaces_ : this->pointTgtFaces_;

    // Determine whether the insertion star is closed or not
    const bool isClosedStar =
        isTri || findIndex(edgeTris_[insertionEdgeis[0]], -1) == -1;

    // Loop the points
    forAll(pointis, pointii)
    {
        const label pointi = pointis[pointii];

        // Initial star tri or edge. If a tri, this means searching the
        // candidates for the one that best contains the point. If an edge then
        // just use the next insertion edge.
        label initialTriOrEdgei = -1;
        if (isTri)
        {
            label insertionCandidateTrii = -1;
            scalar minDistSqr = vGreat;
            forAll(candidateTriTris_, candidateTrii)
            {
                const label trii = candidateTriTris_[candidateTrii];
                if (trii == -1) continue;
                const scalar distSqr = circumDistSqr(trii, pointi);
                if (distSqr < minDistSqr)
                {
                    insertionCandidateTrii = candidateTrii;
                    minDistSqr = distSqr;
                }
                if (minDistSqr < 0)
                {
                    break;
                }
            }
            initialTriOrEdgei = candidateTriTris_[insertionCandidateTrii];
        }
        else
        {
            initialTriOrEdgei = insertionEdgeis[pointii];
        }

        // Populate the star
        star::context starContext = star_.populate
        (
            initialTriOrEdgei,
            isTri,
            [&](const label edgei, const label trii)
            {
                // If this edge is part of a patch edge, or if it connects to a
                // triangle in the other patch, or if it is fixed, then it
                // should not be removed, so it should not be crossed
                if
                (
                    edgePatchEdges(edgei) != labelPair(-1, -1)
                 || (isSrc && triTgtFace_[trii] != -1)
                 || (!isSrc && triSrcFace_[trii] != -1)
                 || findIndex(fixedEdgeis, edgei) != -1
                )
                {
                    return false;
                }

                // If the destination tri borders any tris that are already in
                // the star (other than the tri from which we walked) then do
                // not cross into it. Adding this tri would leave a point
                // within the star, which is not allowed.
                const label triEdgei = findIndex(triEdges_[trii], edgei);
                const label edgej0 = triEdges_[trii][(triEdgei + 1) % 3];
                const label edgej1 = triEdges_[trii][(triEdgei + 2) % 3];
                const label trij0 =
                    edgeTris_[edgej0][edgeTris_[edgej0][0] == trii];
                const label trij1 =
                    edgeTris_[edgej1][edgeTris_[edgej1][0] == trii];
                if
                (
                    (trij0 != -1 && star_.faceStarFaces()[trij0] != -1)
                 || (trij1 != -1 && star_.faceStarFaces()[trij1] != -1)
                 || (!isTri && findIndex(insertionEdgeis, edgej0) != -1)
                 || (!isTri && findIndex(insertionEdgeis, edgej1) != -1)
                )
                {
                    return false;
                }

                // Otherwise, cross into the triangle if it's circumcircle
                // contains the insertion point
                return circumDistSqr(trii, pointi) < 0;
            },
            triEdges_,
            edgeTris_
        );

        // Get the end points of the current insertion edge
        const edge insertionEdge =
            isTri ? edge(-1, -1) : edgePoints(insertionEdgeis[pointii]);

        // Walk around the star polygon disconnecting the old triangles and
        // adding in the new ones
        label newTrii0 = -1, newTrii00 = -1;
        forAllStarEdges(star_, i, starEdgei, edgei)
        {
            const bool isFirst = i == 0;
            const bool isLast = i == star_.starEdgeEdges().size() - 1;

            const bool edgeSidei =
                edgeTris_[edgei][0] == -1
             || star_.faceStarFaces()[edgeTris_[edgei][0]] == -1;

            const label trii = edgeTris_[edgei][edgeSidei];
            const label triEdgei = findIndex(triEdges_[trii], edgei);

            // Disconnect the existing triangle from the star
            const label tempEdgei = newEdgei();
            triEdges_[trii][triEdgei] = tempEdgei;
            edgeTris_[tempEdgei][0] = trii;
            edgeTris_[edgei][edgeSidei] = -1;

            // Add the new triangle
            const label newTrii =
                addTri
                (
                    {
                        pointi,
                        triPoint(trii, triEdgei),
                        triPoint(trii, (triEdgei + 1) % 3),
                    },
                    {
                        !isFirst ? triEdges_[newTrii0][2] : -1,
                        edgei,
                        isClosedStar && isLast ? triEdges_[newTrii00][0] : -1
                    },
                    isSrc ? triSrcFace_[trii] : triTgtFace_[trii],
                    isSrc
                );
            newTrii0 = newTrii;
            newTrii00 = isFirst ? newTrii0 : newTrii00;

            // Add the new triangle into the list of candidates
            if (isTri)
            {
                candidateTriTris_.append(newTrii);
                triCandidateTris_[newTrii] = candidateTriTris_.size() - 1;
            }

            // Get the previous insertion edge
            if (triPoint(newTrii, 2) == insertionEdge[0])
            {
                insertionEdgeis[pointii] = triEdges_[newTrii][2];
            }

            // Get the next insertion edge
            if (!isTri && triPoint(newTrii, 1) == insertionEdge[1])
            {
                insertionEdgeis[pointii + 1] = triEdges_[newTrii][0];
            }
        }

        // Create the edge or face association for the new point
        if (patchEdgei != -1)
        {
            patchEdgePoints[patchEdgei].append(-1);
            label patchEdgePointi = patchEdgePoints[patchEdgei].size() - 1;
            for (; patchEdgePointi >= 0; -- patchEdgePointi)
            {
                label& pointi = patchEdgePoints[patchEdgei][patchEdgePointi];
                pointi = patchEdgePoints[patchEdgei][patchEdgePointi - 1];
                if
                (
                    pointPoints_[pointi] == insertionEdge[0]
                 || pointPoints_[pointi] == insertionEdge[1]
                )
                {
                    break;
                }
            }
            -- patchEdgePointi;
            patchEdgePoints[patchEdgei][patchEdgePointi] = pointi;

            pointPatchEdges[pointi] = patchEdgei;
        }
        else
        {
            pointPatchFaces[pointi] = patchFacei;
        }

        // Check
        if (patchEdgei != -1)
        {
            checkPatchEdge(patchEdgei, isSrc);
        }
        else
        {
            checkPatchFace(patchFacei, isSrc);
        }

        // If the insertion edge is not the same orientation as the original
        // insertion edge then reverse it
        if (!isTri)
        {
            if (edgePoints(insertionEdgeis[pointii + 1])[1] != insertionEdge[1])
            {
                Swap
                (
                    edgeTris_[insertionEdgeis[pointii + 1]][0],
                    edgeTris_[insertionEdgeis[pointii + 1]][1]
                );
            }
        }

        // Remove the star triangles
        forAllStarFaces(star_, starTrii, trii)
        {
            removeTri(trii);
        }

        // Check
        if (patchEdgei != -1)
        {
            checkPatchEdge(patchEdgei, isSrc);
        }
        else
        {
            checkPatchFace(patchFacei, isSrc);
        }
    }

    // Clear the candidate tris
    if (isTri)
    {
        forAll(candidateTriTris_, candidateTrii)
        {
            const label trii = candidateTriTris_[candidateTrii];
            if (trii != -1)
            {
                triCandidateTris_[trii] = -1;
            }
        }
        candidateTriTris_.clear();
    }
}


template<class SrcPatchType, class TgtPatchType>
bool Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::pointCanIntersect
(
    const label pointi
) const
{
    return
        (
            this->pointSrcPoints_[pointi] != -1
         && this->pointTgtPoints_[pointi] == -1
         && this->pointTgtEdges_[pointi] == -1
         && this->pointTgtFaces_[pointi] == -1
        )
     || (
            this->pointTgtPoints_[pointi] != -1
         && this->pointSrcPoints_[pointi] == -1
         && this->pointSrcEdges_[pointi] == -1
         && this->pointSrcFaces_[pointi] == -1
        );
}


template<class SrcPatchType, class TgtPatchType>
bool Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::edgeCanIntersect
(
    const label edgei
) const
{
    return
        (
            edgePatchEdge(edgei, true) != -1
         && edgePatchEdge(edgei, false) == -1
         && (
                edgeTris_[edgei][0] == -1
             || triTgtFace_[edgeTris_[edgei][0]] == -1
            )
         && (
                edgeTris_[edgei][1] == -1
             || triTgtFace_[edgeTris_[edgei][1]] == -1
            )
        )
     || (
            edgePatchEdge(edgei, false) != -1
         && edgePatchEdge(edgei, true) == -1
         && (
                edgeTris_[edgei][0] == -1
             || triSrcFace_[edgeTris_[edgei][0]] == -1
            )
         && (
                edgeTris_[edgei][1] == -1
             || triSrcFace_[edgeTris_[edgei][1]] == -1
            )
        );
}


template<class SrcPatchType, class TgtPatchType>
void
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::snapPatchFaceTris
(
    const label srcFacei,
    const label tgtFacei,
    const scalar snapTol
)
{
    static const FixedList<label, 3> triNoSnap({-1, -1, -1});

    const vector tgtNormal =
        triPointRef(tgtPoints_, patchFacePoints(tgtFacei, false)).normal();

    // Determine what source points snap to
    FixedList<barycentric2D, 3> srcFaceTgtTs(barycentric2D::uniform(-vGreat));
    FixedList<label, 3> srcFaceSnapTgtFacePoint(triNoSnap);
    FixedList<label, 3> srcFaceSnapTgtFaceEdge(triNoSnap);
    forAll(srcFaceSnapTgtFacePoint, srcFacePointi)
    {
        const label srcPointi = patchFacePoint(srcFacei, srcFacePointi, true);

        if (!pointCanIntersect(srcPointi)) continue;

        if ((srcPointNormals_[srcPointi] & tgtNormal) < 0)
        {
            srcFaceTgtTs[srcFacePointi] =
                triIntersect::srcPointTgtTriIntersection
                (
                    srcPoints_[srcPointi],
                    srcPointNormals_[srcPointi],
                    patchFacePointValues(tgtFacei, false, tgtPoints_)
                );

            for (label tgtFacePointi = 0; tgtFacePointi < 3; ++ tgtFacePointi)
            {
                const label tgtPointi =
                    patchFacePoint(tgtFacei, tgtFacePointi, false);

                const label tgtFacePointi0 = (tgtFacePointi + 2) % 3;
                const label tgtFacePointi1 = (tgtFacePointi + 1) % 3;

                if
                (
                    pointCanIntersect(tgtPointi)
                 && mag(srcFaceTgtTs[srcFacePointi][tgtFacePointi0]) < snapTol
                 && mag(srcFaceTgtTs[srcFacePointi][tgtFacePointi1]) < snapTol
                )
                {
                    srcFaceSnapTgtFacePoint[srcFacePointi] =
                        srcFaceSnapTgtFacePoint[srcFacePointi] == -1
                      ? tgtFacePointi
                      : -2;

                    srcFaceTgtTs[srcFacePointi][tgtFacePointi] = 1;
                    srcFaceTgtTs[srcFacePointi][tgtFacePointi0] = 0;
                    srcFaceTgtTs[srcFacePointi][tgtFacePointi1] = 0;
                }
            }

            srcFaceSnapTgtFacePoint[srcFacePointi] =
                max(srcFaceSnapTgtFacePoint[srcFacePointi], -1);

            if (srcFaceSnapTgtFacePoint[srcFacePointi] != -1) continue;

            for (label tgtFaceEdgei = 0; tgtFaceEdgei < 3; ++ tgtFaceEdgei)
            {
                const label tgtFacePointi0 = tgtFaceEdgei;
                const label tgtFacePointi1 = (tgtFaceEdgei + 1) % 3;
                const label tgtFacePointiOpp = (tgtFaceEdgei + 2) % 3;

                if
                (
                    srcFaceTgtTs[srcFacePointi][tgtFacePointi0] >= snapTol
                 && srcFaceTgtTs[srcFacePointi][tgtFacePointi1] >= snapTol
                 && mag(srcFaceTgtTs[srcFacePointi][tgtFacePointiOpp]) < snapTol
                )
                {
                    srcFaceSnapTgtFaceEdge[srcFacePointi] =
                        srcFaceSnapTgtFaceEdge[srcFacePointi] == -1
                      ? tgtFaceEdgei
                      : -2;

                    const scalar eps =
                        srcFaceTgtTs[srcFacePointi][tgtFacePointiOpp];
                    srcFaceTgtTs[srcFacePointi][tgtFacePointiOpp] = 0;
                    srcFaceTgtTs[srcFacePointi][tgtFacePointi0] /= 1 - eps;
                    srcFaceTgtTs[srcFacePointi][tgtFacePointi1] /= 1 - eps;
                }
            }

            srcFaceSnapTgtFaceEdge[srcFacePointi] =
                max(srcFaceSnapTgtFaceEdge[srcFacePointi], -1);
        }
    }

    // Determine what target points snap to
    FixedList<barycentric2D, 3> tgtFaceSrcTs(barycentric2D::uniform(-vGreat));
    FixedList<label, 3> tgtFaceSnapSrcFacePoint(triNoSnap);
    FixedList<label, 3> tgtFaceSnapSrcFaceEdge(triNoSnap);
    forAll(tgtFaceSnapSrcFacePoint, tgtFacePointi)
    {
        const label tgtPointi = patchFacePoint(tgtFacei, tgtFacePointi, false);

        if (!pointCanIntersect(tgtPointi)) continue;

        tgtFaceSrcTs[tgtFacePointi] =
            triIntersect::srcTriTgtPointIntersection
            (
                patchFacePointValues(srcFacei, true, srcPoints_),
                patchFacePointValues(srcFacei, true, srcPointNormals_),
                tgtPoints_[tgtPointi]
            );

        const vector srcNormal =
            triIntersect::srcTriInterpolate
            (
                tgtFaceSrcTs[tgtFacePointi],
                patchFacePointValues(srcFacei, true, srcPointNormals_)
            );

        if ((srcNormal & tgtNormal) < 0)
        {
            for (label srcFacePointi = 0; srcFacePointi < 3; ++ srcFacePointi)
            {
                const label srcPointi =
                    patchFacePoint(srcFacei, srcFacePointi, true);

                const label srcFacePointi0 = (srcFacePointi + 2) % 3;
                const label srcFacePointi1 = (srcFacePointi + 1) % 3;

                if
                (
                    pointCanIntersect(srcPointi)
                 && mag(tgtFaceSrcTs[tgtFacePointi][srcFacePointi0]) < snapTol
                 && mag(tgtFaceSrcTs[tgtFacePointi][srcFacePointi1]) < snapTol
                )
                {
                    tgtFaceSnapSrcFacePoint[tgtFacePointi] =
                        tgtFaceSnapSrcFacePoint[tgtFacePointi] == -1
                      ? srcFacePointi
                      : -2;

                    tgtFaceSrcTs[tgtFacePointi][srcFacePointi] = 1;
                    tgtFaceSrcTs[tgtFacePointi][srcFacePointi0] = 0;
                    tgtFaceSrcTs[tgtFacePointi][srcFacePointi1] = 0;
                }
            }

            tgtFaceSnapSrcFacePoint[tgtFacePointi] =
                max(tgtFaceSnapSrcFacePoint[tgtFacePointi], -1);

            if (tgtFaceSnapSrcFacePoint[tgtFacePointi] != -1) continue;

            for (label srcFaceEdgei = 0; srcFaceEdgei < 3; ++ srcFaceEdgei)
            {
                const label srcFacePointi0 = srcFaceEdgei;
                const label srcFacePointi1 = (srcFaceEdgei + 1) % 3;
                const label srcFacePointiOpp = (srcFaceEdgei + 2) % 3;

                if
                (
                    tgtFaceSrcTs[tgtFacePointi][srcFacePointi0] >= snapTol
                 && tgtFaceSrcTs[tgtFacePointi][srcFacePointi1] >= snapTol
                 && mag(tgtFaceSrcTs[tgtFacePointi][srcFacePointiOpp]) < snapTol
                )
                {
                    tgtFaceSnapSrcFaceEdge[tgtFacePointi] =
                        tgtFaceSnapSrcFaceEdge[tgtFacePointi] == -1
                      ? srcFaceEdgei
                      : -2;

                    const scalar eps =
                        tgtFaceSrcTs[tgtFacePointi][srcFacePointiOpp];
                    tgtFaceSrcTs[tgtFacePointi][srcFacePointiOpp] = 0;
                    tgtFaceSrcTs[tgtFacePointi][srcFacePointi0] /= 1 - eps;
                    tgtFaceSrcTs[tgtFacePointi][srcFacePointi1] /= 1 - eps;
                }
            }

            tgtFaceSnapSrcFaceEdge[tgtFacePointi] =
                max(tgtFaceSnapSrcFaceEdge[tgtFacePointi], -1);
        }
    }

    // Return if there is nothing to do
    if
    (
        srcFaceSnapTgtFacePoint == triNoSnap
     && srcFaceSnapTgtFaceEdge == triNoSnap
     && tgtFaceSnapSrcFacePoint == triNoSnap
     && tgtFaceSnapSrcFaceEdge == triNoSnap
    )
    {
        return;
    }

    // Make point-point snaps symmetric
    forAll(srcFaceSnapTgtFacePoint, srcFacePointi)
    {
        const label tgtFacePointi = srcFaceSnapTgtFacePoint[srcFacePointi];

        if (tgtFacePointi != -1)
        {
            tgtFaceSnapSrcFacePoint[tgtFacePointi] = srcFacePointi;
        }
    }
    forAll(tgtFaceSnapSrcFacePoint, tgtFacePointi)
    {
        const label srcFacePointi = tgtFaceSnapSrcFacePoint[tgtFacePointi];

        if (srcFacePointi != -1)
        {
            srcFaceSnapTgtFacePoint[srcFacePointi] = tgtFacePointi;
        }
    }

    // Do point-point motion
    forAll(srcFaceSnapTgtFacePoint, srcFacePointi)
    {
        const label tgtFacePointi = srcFaceSnapTgtFacePoint[srcFacePointi];

        if (tgtFacePointi == -1) continue;

        const label srcPointi = patchFacePoint(srcFacei, srcFacePointi, true);
        const label tgtPointi = patchFacePoint(tgtFacei, tgtFacePointi, false);

        if (!pointCanIntersect(srcPointi)) continue;
        if (!pointCanIntersect(tgtPointi)) continue;

        srcPoints_[tgtPointi] = srcPoints_[srcPointi];
        srcPointNormals_[tgtPointi] = srcPointNormals_[srcPointi];
        tgtPoints_[srcPointi] = tgtPoints_[tgtPointi];

        const point& srcP = srcPoints_[srcPointi];
        const vector& srcN = srcPointNormals_[srcPointi];
        const point& tgtP = tgtPoints_[tgtPointi];
        const vector d = (tensor::I - sqr(srcN)) & (tgtP - srcP);
        srcPoints_[srcPointi] += d/2;
        tgtPoints_[tgtPointi] -= d/2;
    }

    // Do source-point-target-edge motion
    forAll(srcFaceSnapTgtFaceEdge, srcFacePointi)
    {
        const label tgtFaceEdgei = srcFaceSnapTgtFaceEdge[srcFacePointi];

        if (tgtFaceEdgei == -1) continue;

        const label srcPointi = patchFacePoint(srcFacei, srcFacePointi, true);

        if (!pointCanIntersect(srcPointi)) continue;

        tgtPoints_[srcPointi] =
            triIntersect::tgtTriInterpolate
            (
                srcFaceTgtTs[srcFacePointi],
                patchFacePointValues(tgtFacei, false, tgtPoints_)
            );

        const point& srcP = srcPoints_[srcPointi];
        const vector& srcN = srcPointNormals_[srcPointi];
        const point& tgtP = tgtPoints_[srcPointi];
        const vector d = (tensor::I - sqr(srcN)) & (tgtP - srcP);
        srcPoints_[srcPointi] += d;
    }

    // Do target-point-source-edge motion
    forAll(tgtFaceSnapSrcFaceEdge, tgtFacePointi)
    {
        const label srcFaceEdgei = tgtFaceSnapSrcFaceEdge[tgtFacePointi];

        if (srcFaceEdgei == -1) continue;

        const label tgtPointi = patchFacePoint(tgtFacei, tgtFacePointi, false);

        if (!pointCanIntersect(tgtPointi)) continue;

        srcPoints_[tgtPointi] =
            triIntersect::srcTriInterpolate
            (
                tgtFaceSrcTs[tgtFacePointi],
                patchFacePointValues(srcFacei, true, srcPoints_)
            );
        srcPointNormals_[tgtPointi] =
            triIntersect::srcTriInterpolate
            (
                tgtFaceSrcTs[tgtFacePointi],
                patchFacePointValues(srcFacei, true, srcPointNormals_)
            );

        const point& srcP = srcPoints_[tgtPointi];
        const vector& srcN = srcPointNormals_[tgtPointi];
        const point& tgtP = tgtPoints_[tgtPointi];
        const vector d = (tensor::I - sqr(srcN)) & (tgtP - srcP);
        tgtPoints_[tgtPointi] -= d;
    }

    // Do point-point snapping
    forAll(srcFaceSnapTgtFacePoint, srcFacePointi)
    {
        const label tgtFacePointi = srcFaceSnapTgtFacePoint[srcFacePointi];

        if (tgtFacePointi == -1) continue;

        const label srcPointi = patchFacePoint(srcFacei, srcFacePointi, true);
        const label tgtPointi = patchFacePoint(tgtFacei, tgtFacePointi, false);

        if (!pointCanIntersect(srcPointi)) continue;
        if (!pointCanIntersect(tgtPointi)) continue;

        srcPoints_[tgtPointi] = srcPoints_[srcPointi];
        srcPointNormals_[tgtPointi] = srcPointNormals_[srcPointi];
        tgtPoints_[srcPointi] = tgtPoints_[tgtPointi];

        pointPoints_[tgtPointi] = srcPointi;

        this->pointTgtEdges_[srcPointi] = this->pointTgtEdges_[tgtPointi];
        this->pointSrcEdges_[tgtPointi] = -1;
        this->pointTgtEdges_[tgtPointi] = -1;

        this->pointTgtPoints_[srcPointi] = this->pointTgtPoints_[tgtPointi];
        this->tgtPointPoints_[this->pointTgtPoints_[tgtPointi]] = srcPointi;
        this->pointSrcPoints_[tgtPointi] = -1;
        this->pointTgtPoints_[tgtPointi] = -1;
    }

    // Check
    checkPatchFace(srcFacei, true);
    checkPatchFace(tgtFacei, false);

    // Insertion workspace
    labelList insertPointis(1), insertEdgeis(2, label(-1)), fixedEdgeis(0);

    // Do source-point-target-edge snapping
    forAll(srcFaceSnapTgtFaceEdge, srcFacePointi)
    {
        const label tgtFaceEdgei = srcFaceSnapTgtFaceEdge[srcFacePointi];

        if (tgtFaceEdgei == -1) continue;

        const label srcPointi = patchFacePoint(srcFacei, srcFacePointi, true);

        if (!pointCanIntersect(srcPointi)) continue;

        const label tgtPatchEdgei =
            patchFacePatchEdge(tgtFacei, tgtFaceEdgei, false);

        // Loop all target face tris and find the best edge into which to
        // insert the source point
        label insertTgtEdgei = -1;
        scalar insertDistSqr = vGreat;
        forAll(tgtFaceTris_[tgtFacei], tgtFaceTrii)
        {
            const label tgtTrii = tgtFaceTris_[tgtFacei][tgtFaceTrii];

            forAll(triEdges_[tgtTrii], tgtTriEdgei)
            {
                const label tgtEdgei = triEdges_[tgtTrii][tgtTriEdgei];

                if
                (
                    edgeCanIntersect(tgtEdgei)
                 && tgtPatchEdgei == edgePatchEdge(tgtEdgei, false)
                )
                {
                    const scalar distSqr = circumDistSqr(tgtTrii, srcPointi);

                    if (distSqr < insertDistSqr)
                    {
                        insertDistSqr = distSqr;
                        insertTgtEdgei = tgtEdgei;
                    }
                }
            }
        }

        if (insertTgtEdgei != -1)
        {
            insertPointis[0] = srcPointi;
            insertPoints
            (
                insertTgtEdgei,
                false,
                insertPointis,
                insertEdgeis,
                fixedEdgeis
            );
        }
    }

    // Do target-point-source-edge snapping
    forAll(tgtFaceSnapSrcFaceEdge, tgtFacePointi)
    {
        const label srcFaceEdgei = tgtFaceSnapSrcFaceEdge[tgtFacePointi];

        if (srcFaceEdgei == -1) continue;

        const label tgtPointi = patchFacePoint(tgtFacei, tgtFacePointi, false);

        if (!pointCanIntersect(tgtPointi)) continue;

        const label srcPatchEdgei =
            patchFacePatchEdge(srcFacei, srcFaceEdgei, true);

        srcPoints_[tgtPointi] =
            triIntersect::srcTriInterpolate
            (
                tgtFaceSrcTs[tgtFacePointi],
                patchFacePointValues(srcFacei, true, srcPoints_)
            );
        srcPointNormals_[tgtPointi] =
            triIntersect::srcTriInterpolate
            (
                tgtFaceSrcTs[tgtFacePointi],
                patchFacePointValues(srcFacei, true, srcPointNormals_)
            );

        // Loop all target face tris and find the best edge into which to
        // insert the source point
        label insertSrcEdgei = -1;
        scalar insertDistSqr = vGreat;
        forAll(srcFaceTris_[srcFacei], srcFaceTrii)
        {
            const label srcTrii = srcFaceTris_[srcFacei][srcFaceTrii];

            forAll(triEdges_[srcTrii], srcTriEdgei)
            {
                const label srcEdgei = triEdges_[srcTrii][srcTriEdgei];

                if
                (
                    edgeCanIntersect(srcEdgei)
                 && srcPatchEdgei == edgePatchEdge(srcEdgei, true)
                )
                {
                    const scalar distSqr = circumDistSqr(srcTrii, tgtPointi);

                    if (distSqr < insertDistSqr)
                    {
                        insertDistSqr = distSqr;
                        insertSrcEdgei = srcEdgei;
                    }
                }
            }
        }

        if (insertSrcEdgei != -1)
        {
            insertPointis[0] = tgtPointi;
            insertPoints
            (
                insertSrcEdgei,
                false,
                insertPointis,
                insertEdgeis,
                fixedEdgeis
            );
        }
    }

    // Check
    checkPatchFace(srcFacei, true);
    checkPatchFace(tgtFacei, false);
}


template<class SrcPatchType, class TgtPatchType>
bool Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::intersectTris
(
    const label srcTrii,
    const label tgtTrii
)
{
    // Get the source and target faces that these tris are associated with
    const label srcFacei = triSrcFace_[srcTrii];
    const label tgtFacei = triTgtFace_[tgtTrii];
    if (srcFacei == -1 || tgtFacei == -1)
    {
        FatalErrorInFunction
            << "Tri-intersections must be between a tri associated with the "
            << "source patch and one associated with the target patch"
            << exit(FatalError);
    }

    // Do intersection
    DynamicList<point> ictSrcPoints;
    DynamicList<vector> ictSrcPointNormals;
    DynamicList<point> ictTgtPoints;
    DynamicList<triIntersect::location> ictPointLocations;
    triIntersect::intersectTris
    (
        triPointValues(srcTrii, srcPoints_),
        triPointValues(srcTrii, srcPointNormals_),
        triOwns(srcTrii),
        triOtherTriPoints(srcTrii, tgtTrii),
        triPointValues(tgtTrii, tgtPoints_),
        triOwns(tgtTrii),
        triOtherTriPoints(tgtTrii, srcTrii),
        ictSrcPoints,
        ictSrcPointNormals,
        ictTgtPoints,
        ictPointLocations,
        this->debug > 3,
        this->debug > 4
      ? "srcTrii##tgtTrii=" + name(srcTrii) + name(tgtTrii)
      : ""
    );

    // Return false if no intersection
    if (!ictPointLocations.size())
    {
        return false;
    }

    // Loop over the intersection points looking for elements that are already
    // intersected with the opposing surface. If any are found, then this
    // intersection is considered invalid.
    forAll(ictPointLocations, ictPointi)
    {
        const triIntersect::location& l = ictPointLocations[ictPointi];

        if (l.isSrcPoint() && !l.isTgtPoint())
        {
            const label pointi = triPoint(srcTrii, l.srcPointi());

            if
            (
                this->pointTgtPoints_[pointi] != -1
             || this->pointTgtEdges_[pointi] != -1
             || this->pointTgtFaces_[pointi] != -1
            )
            {
                return false;
            }
        }
        if (!l.isSrcPoint() && l.isTgtPoint())
        {
            const label pointi = triPoint(tgtTrii, l.tgtPointi());

            if
            (
                this->pointSrcPoints_[pointi] != -1
             || this->pointSrcEdges_[pointi] != -1
             || this->pointSrcFaces_[pointi] != -1
            )
            {
                return false;
            }
        }
        if (l.isIntersection())
        {
            const label srcEdgei = triEdges_[srcTrii][l.srcEdgei()];
            const label tgtEdgei = triEdges_[tgtTrii][l.tgtEdgei()];

            const labelPair srcPatchEdges = edgePatchEdges(srcEdgei);
            const labelPair tgtPatchEdges = edgePatchEdges(tgtEdgei);

            if
            (
                (srcPatchEdges[0] != -1 && tgtPatchEdges[1] != -1)
             && (srcPatchEdges[1] != -1 || tgtPatchEdges[0] != -1)
            )
            {
                return false;
            }

            const label srcTrij =
                edgeTris_[srcEdgei][edgeTris_[srcEdgei][0] == srcTrii];
            const label tgtTrij =
                edgeTris_[tgtEdgei][edgeTris_[tgtEdgei][0] == tgtTrii];

            if
            (
                (srcTrij != -1 && triTgtFace_[srcTrij] != -1)
             || (tgtTrij != -1 && triSrcFace_[tgtTrij] != -1)
            )
            {
                return false;
            }
        }
    }

    // Flags controlling what operations to do
    bool insertPointsIntoTri = false;
    bool insertPointsIntoEdges = false;

    // Loop over the intersection points creating a map from intersection point
    // indices to point indices. Overwrite any intersected source or target
    // points, and create points generated by any new intersections. Point-tri
    // and point-edge associations are generated later during insertion.
    labelList ictPointiToPointi(ictPointLocations.size(), -1);
    forAll(ictPointLocations, ictPointi)
    {
        const triIntersect::location& l = ictPointLocations[ictPointi];

        if (l.isSrcPoint() && !l.isTgtPoint())
        {
            const label pointi = triPoint(srcTrii, l.srcPointi());

            ictPointiToPointi[ictPointi] = pointi;

            // Store the projected location on the target side
            tgtPoints_[pointi] = ictTgtPoints[ictPointi];

            insertPointsIntoTri = true;
        }
        if (!l.isSrcPoint() && l.isTgtPoint())
        {
            const label pointi = triPoint(tgtTrii, l.tgtPointi());

            ictPointiToPointi[ictPointi] = pointi;

            // Store the projected location on the source side
            srcPoints_[pointi] = ictSrcPoints[ictPointi];
            srcPointNormals_[pointi] = ictSrcPointNormals[ictPointi];

            insertPointsIntoTri = true;
        }
        if (l.isIntersection())
        {
            const label srcPatchEdgei =
                edgePatchEdge(triEdges_[srcTrii][l.srcEdgei()], true);
            const label tgtPatchEdgei =
                edgePatchEdge(triEdges_[tgtTrii][l.tgtEdgei()], false);

            if (srcPatchEdgei != -1 && tgtPatchEdgei != -1)
            {
                ictPointiToPointi[ictPointi] = srcPoints_.size();

                const label pointi = newPointi();

                srcPoints_[pointi] = ictSrcPoints[ictPointi];
                srcPointNormals_[pointi] = ictSrcPointNormals[ictPointi];
                tgtPoints_[pointi] = ictTgtPoints[ictPointi];

                insertPointsIntoEdges = true;
            }
        }
    }

    // Store the source and target tri edges and edge-ownerships
    const FixedList<label, 3> srcTriEdges = triEdges_[srcTrii];
    const FixedList<bool, 3> srcTriOwns = triOwns(srcTrii);
    const FixedList<label, 3> tgtTriEdges = triEdges_[tgtTrii];
    const FixedList<bool, 3> tgtTriOwns = triOwns(tgtTrii);

    // Edges which are to be fixed during insertion. These constraints are not
    // needed. The only edges that get inserted into are those that associate
    // with patch edges, and these are constrained anyway. Constraining all
    // edges of the tris leads to poor triangle qualities, as it limits the
    // that re-triangulation insertPoints can do. For testing, however, this
    // can be quite useful, hence the option to make the system over-constrain
    // the problem.
    DynamicList<label> fixedSrcEdgeis, fixedTgtEdgeis;
    #ifdef OVERCONSTRAIN
    forAll(srcTriEdges, srcTriEdgei)
    {
        fixedSrcEdgeis.append(srcTriEdges[srcTriEdgei]);
    }
    forAll(tgtTriEdges, tgtTriEdgei)
    {
        fixedTgtEdgeis.append(tgtTriEdges[tgtTriEdgei]);
    }
    #endif

    // Insert points into the triangles
    if (insertPointsIntoTri)
    {
        DynamicList<label> insertSrcPointis(3), insertTgtPointis(3);
        DynamicList<label> insertEdgeis;
        forAll(ictPointLocations, ictPointi)
        {
            const label pointi = ictPointiToPointi[ictPointi];

            if (pointi == -1) continue;

            const triIntersect::location& l = ictPointLocations[ictPointi];

            if (l.isSrcPoint() && !l.isTgtPoint())
            {
                insertSrcPointis.append(pointi);
            }
            if (!l.isSrcPoint() && l.isTgtPoint())
            {
                insertTgtPointis.append(pointi);
            }
        }

        if (insertSrcPointis.size())
        {
            insertPoints
            (
                tgtTrii,
                true,
                insertSrcPointis,
                insertEdgeis,
                fixedTgtEdgeis
            );
        }

        if (insertTgtPointis.size())
        {
            insertPoints
            (
                srcTrii,
                true,
                insertTgtPointis,
                insertEdgeis,
                fixedSrcEdgeis
            );
        }

        // Check
        checkPatchFace(srcFacei, true);
        checkPatchFace(tgtFacei, false);
    }

    // Insert points into the edges
    if (insertPointsIntoEdges)
    {
        DynamicList<label> insertPointis(2), insertEdgeis(3);
        forAll(ictPointLocations, ictPointi0)
        {
            const label ictPointi1 = ictPointLocations.fcIndex(ictPointi0);

            const label pointi0 = ictPointiToPointi[ictPointi0];
            const label pointi1 = ictPointiToPointi[ictPointi1];

            const triIntersect::location& l0 = ictPointLocations[ictPointi0];
            const triIntersect::location& l1 = ictPointLocations[ictPointi1];

            insertPointis.clear();
            insertEdgeis.clear();
            insertEdgeis.append(-1);

            if
            (
                (
                    l0.isIntersection()
                 && l1.isSrcPoint()
                 && l0.srcEdgei() != (l1.srcPointi() + 1) % 3
                )
             || (
                    l0.isSrcPoint()
                 && l1.isIntersection()
                 && (l0.srcPointi() + 1) % 3 != l1.srcEdgei()
                )
             || (
                    l0.isIntersection()
                 && l1.isIntersection()
                 && l0.srcEdgei() == l1.srcEdgei()
                )
            )
            {
                const label srcTriEdgei =
                    l0.isIntersection() ? l0.srcEdgei() : l1.srcEdgei();

                if (!l0.isSrcPoint() && pointi0 != -1)
                {
                    insertPointis.append(pointi0);
                    insertEdgeis.append(-1);
                }
                if (!l1.isSrcPoint() && pointi1 != -1)
                {
                    insertPointis.append(pointi1);
                    insertEdgeis.append(-1);
                }
                if (!srcTriOwns[srcTriEdgei])
                {
                    inplaceReverseList(insertPointis);
                }

                if (insertPointis.size())
                {
                    insertPoints
                    (
                        srcTriEdges[srcTriEdgei],
                        false,
                        insertPointis,
                        insertEdgeis,
                        fixedSrcEdgeis
                    );
                }

                #ifdef OVERCONSTRAIN
                fixedSrcEdgeis[srcTriEdgei] = insertEdgeis.remove();
                fixedSrcEdgeis.append(insertEdgeis);
                #endif
            }

            if
            (
                (
                    l0.isIntersection()
                 && l1.isTgtPoint()
                 && l0.tgtEdgei() != (l1.tgtPointi() + 1) % 3
                )
             || (
                    l0.isTgtPoint()
                 && l1.isIntersection()
                 && (l0.tgtPointi() + 1) % 3 != l1.tgtEdgei()
                )

             || (
                    l0.isIntersection()
                 && l1.isIntersection()
                 && l0.tgtEdgei() == l1.tgtEdgei()
                )
            )
            {
                const label tgtTriEdgei =
                    l0.isIntersection() ? l0.tgtEdgei() : l1.tgtEdgei();

                if (!l0.isTgtPoint() && pointi0 != -1)
                {
                    insertPointis.append(pointi0);
                    insertEdgeis.append(-1);
                }
                if (!l1.isTgtPoint() && pointi1 != -1)
                {
                    insertPointis.append(pointi1);
                    insertEdgeis.append(-1);
                }
                if (tgtTriOwns[tgtTriEdgei])
                {
                    inplaceReverseList(insertPointis);
                }

                if (insertPointis.size())
                {
                    insertPoints
                    (
                        tgtTriEdges[tgtTriEdgei],
                        false,
                        insertPointis,
                        insertEdgeis,
                        fixedTgtEdgeis
                    );
                }

                #ifdef OVERCONSTRAIN
                fixedTgtEdgeis[tgtTriEdgei] = insertEdgeis.remove();
                fixedTgtEdgeis.append(insertEdgeis);
                #endif
            }
        }

        // Check
        checkPatchFace(srcFacei, true);
        checkPatchFace(tgtFacei, false);
    }

    return true;
}


template<class SrcPatchType, class TgtPatchType>
void
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::intersectPatchFaceTris
(
    const label srcFacei,
    const label tgtFacei
)
{
    // Get the source face tris and make sure none are marked
    const DynamicList<label>& srcFaceTris = srcFaceTris_[srcFacei];
    forAll(srcFaceTris, srcFaceTrii)
    {
        const label srcTrii = srcFaceTris[srcFaceTrii];

        if (triMarkedTris_[srcTrii] != -1)
        {
            markedTriTris_[triMarkedTris_[srcTrii]] = -1;
            triMarkedTris_[srcTrii] = -1;
        }
    }

    // Loop the source face tris
    label srcFaceTrii = 0;
    while (srcFaceTrii < srcFaceTris.size())
    {
        // Mark this source tri
        const label srcTrii = srcFaceTris[srcFaceTrii];
        markedTriTris_.append(srcTrii);
        triMarkedTris_[srcTrii] = markedTriTris_.size() - 1;

        // Get the target face tris and make sure none are marked
        const DynamicList<label>& tgtFaceTris = tgtFaceTris_[tgtFacei];
        forAll(tgtFaceTris, tgtFaceTrii)
        {
            const label tgtTrii = tgtFaceTris[tgtFaceTrii];

            if (triMarkedTris_[tgtTrii] != -1)
            {
                markedTriTris_[triMarkedTris_[tgtTrii]] = -1;
                triMarkedTris_[tgtTrii] = -1;
            }
        }

        // Loop the target face tris
        label tgtFaceTrii = 0;
        while (tgtFaceTrii < tgtFaceTris.size())
        {
            // Mark this target tri
            const label tgtTrii = tgtFaceTris[tgtFaceTrii];
            markedTriTris_.append(tgtTrii);
            triMarkedTris_[tgtTrii] = markedTriTris_.size() - 1;

            /*
            // If the target tri has not unmarked then snap
            if (triMarkedTris_[tgtFaceTris[tgtFaceTrii]] != -1)
            {
                snapTris(srcTrii, tgtTrii);
            }

            // If the source tri has unmarked then break
            if (triMarkedTris_[srcFaceTris[srcFaceTrii]] == -1) break;
            */

            // If the target tri has not unmarked then intersect
            if (triMarkedTris_[tgtFaceTris[tgtFaceTrii]] != -1)
            {
                intersectTris(srcTrii, tgtTrii);
            }

            // If the source tri has unmarked then break
            if (triMarkedTris_[srcFaceTris[srcFaceTrii]] == -1) break;

            // Iterate backwards to find the next marked target tri
            while
            (
                tgtFaceTrii >= 0
             && triMarkedTris_[tgtFaceTris[tgtFaceTrii]] == -1
            )
            {
                -- tgtFaceTrii;
            }

            // Increment to the next unmarked target tri
            ++ tgtFaceTrii;
        }

        // Iterate backwards to find the next marked source tri
        while
        (
            srcFaceTrii >= 0
         && triMarkedTris_[srcFaceTris[srcFaceTrii]] == -1
        )
        {
            -- srcFaceTrii;
        }

        // Increment to the next unmarked source tri
        ++ srcFaceTrii;
    }

    // Clear the marked tris
    forAll(markedTriTris_, markedTrii)
    {
        const label trii = markedTriTris_[markedTrii];
        if (trii != -1)
        {
            triMarkedTris_[trii] = -1;
        }
    }
    markedTriTris_.clear();
}


template<class SrcPatchType, class TgtPatchType>
bool
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::conformPatchFaceTris
(
    const label patchFacei,
    const label otherPatchFacei,
    const bool isSrc
)
{
    const labelList& patchFaceTris =
        isSrc ? srcFaceTris_[patchFacei] : tgtFaceTris_[patchFacei];
    const labelList& otherPatchFaceTris =
        isSrc ? tgtFaceTris_[otherPatchFacei] : srcFaceTris_[otherPatchFacei];

    const labelList& triOtherPatchFaces =
        isSrc ? triTgtFace_ : triSrcFace_;

    const pointField& points = isSrc ? srcPoints_ : tgtPoints_;

    const vector& patchFaceN =
        isSrc
      ? this->srcPatch_.faceNormals()[patchFacei]
      : this->tgtPatch_.faceNormals()[patchFacei];

    auto pointIntersectsPatchFace = [&]
    (
        const label pointi,
        const label patchFacei,
        const bool isSrc
    )
    {
        const labelList& pointPatchPoints =
            isSrc ? this->pointSrcPoints_ : this->pointTgtPoints_;
        const labelList& pointPatchEdges =
            isSrc ? this->pointSrcEdges_ : this->pointTgtEdges_;
        const labelList& pointPatchFaces =
            isSrc ? this->pointSrcFaces_ : this->pointTgtFaces_;

        const triFace patchFacePatchPoints
        (
            this->patchFacePatchPoints(patchFacei, isSrc)
        );
        const triFace patchFacePatchEdges
        (
            this->patchFacePatchEdges(patchFacei, isSrc)
        );

        return
            findIndex
            (
                patchFacePatchPoints,
                pointPatchPoints[pointi]
            ) != -1
         || findIndex
            (
                patchFacePatchEdges,
                pointPatchEdges[pointi]
            ) != -1
         || pointPatchFaces[pointi] == patchFacei;
    };

    auto edgeIntersectsPatchFace = [&]
    (
        const label edgei,
        const label patchFacei,
        const bool isSrc
    )
    {
        const edge e = edgePoints(edgei);
        return
            pointIntersectsPatchFace(e[0], patchFacei, isSrc)
         && pointIntersectsPatchFace(e[1], patchFacei, isSrc);
    };

    // Get all the edges that need conforming to
    DynamicList<label> conformEdgeis;
    forAll(otherPatchFaceTris, otherPatchFaceTrii)
    {
        const label trii = otherPatchFaceTris[otherPatchFaceTrii];

        forAll(triEdges_[trii], triEdgei)
        {
            const label edgei = triEdges_[trii][triEdgei];

            const label trij =
                edgeTris_[edgei][edgeTris_[edgei][0] == trii];

            const label otherPatchFacej =
                trij == -1 ? -1 : triOtherPatchFaces[trij];

            if
            (
                otherPatchFacei != otherPatchFacej
             && edgeIntersectsPatchFace(edgei, patchFacei, isSrc)
            )
            {
                conformEdgeis.append(edgei);
            }
        }
    }

    // Remove edges already present in the triangulation by shuffling up
    auto isConformed = [&](const label edgei)
    {
        const edge e = edgePoints(edgei);

        forAll(patchFaceTris, patchFaceTrii)
        {
            const label trii = patchFaceTris[patchFaceTrii];

            forAll(triEdges_[trii], triEdgei)
            {
                const label edgei = triEdges_[trii][triEdgei];

                if (edge::compare(e, edgePoints(edgei)) != 0)
                {
                    return true;
                }
            }
        }

        return false;
    };
    {
        label conformi = 0;
        forAll(conformEdgeis, conformj)
        {
            if (!isConformed(conformEdgeis[conformj]))
            {
                conformEdgeis[conformi] = conformEdgeis[conformj];
                ++ conformi;
            }
        }
        conformEdgeis.resize(conformi);
    }

    // Quit if there is nothing to conform to
    if (conformEdgeis.empty()) return true;

    // Get the next patch face tri index connected to a given point
    auto nextConnectedPatchFaceTri = [&]
    (
        const label pointi,
        const label patchFaceTrii0 = 0
    )
    {
        for
        (
            label patchFaceTrii = patchFaceTrii0;
            patchFaceTrii < patchFaceTris.size();
            ++ patchFaceTrii
        )
        {
            const label trii = patchFaceTris[patchFaceTrii];

            if (findIndex(triPoints(trii), pointi) != -1)
            {
                return patchFaceTrii;
            }
        }

        return label(-1);
    };

    // Track from a point on a triangle towards a given point. Stop at an edge
    // and set the index of that edge and update the local coordinates.
    auto trackToEdge = [&]
    (
        const label trii,
        label& triEdgei,
        barycentric2D& y,
        const label pointi1
    )
    {
        const triPointRef t = triPoints(trii).tri(points);

        const barycentricTensor2D A(t.a(), t.b(), t.c());

        const point p = A & y;
        const vector dp = points[pointi1] - p;

        const vector ab = t.b() - t.a();
        const vector ac = t.c() - t.a();
        const vector bc = t.c() - t.b();
        const scalar detA = (ab ^ ac) & patchFaceN;
        const barycentricTensor2D T
        (
            patchFaceN ^ bc,
            ac ^ patchFaceN,
            patchFaceN ^ ab
        );

        const barycentric2D TDp = dp & T;

        label iH = -1;
        scalar lambdaByDetAH =
            !std::isnormal(detA) || detA < 0 ? vGreat : 1/detA;

        forAll(TDp, i)
        {
            if (TDp[i] < - detA*small)
            {
                const scalar lambdaByDetA = - y[i]/TDp[i];

                if (0 <= lambdaByDetA && lambdaByDetA < lambdaByDetAH)
                {
                    iH = i;
                    lambdaByDetAH = lambdaByDetA;
                }
            }
        }

        y += lambdaByDetAH*TDp;
        forAll(y, i)
        {
            y.replace(i, i == iH ? 0 : max(0, y[i]));
        }
        if (iH == -1)
        {
            y /= cmptSum(y);
        }

        triEdgei = iH == -1 ? -1 : (iH + 1) % 3;
    };

    // Cross the edge. Return true if the edge can be crossed, and false
    // otherwise. Set the new triangle, new edge index, and transform the local
    // coordinates appropriately.
    auto crossEdge = [&](label& trii, label& triEdgei, barycentric2D& y)
    {
        if (triEdgei == -1) return false;

        const label edgei = triEdges_[trii][triEdgei];

        if (edgePatchEdges(edgei) != labelPair(-1, -1)) return false;

        const label trij = edgeTris_[edgei][edgeTris_[edgei][0] == trii];

        if (trij == -1) return false;

        if (triSrcFace_[trii] != triSrcFace_[trij]) return false;
        if (triTgtFace_[trii] != triTgtFace_[trij]) return false;

        const label triEdgej = findIndex(triEdges_[trij], edgei);

        auto inplaceRotate = [](barycentric2D& y, label n)
        {
            n = n % 3;

            if (n == 1)
            {
                scalar temp = y.a();
                y.a() = y.b();
                y.b() = y.c();
                y.c() = temp;
            }

            if (n == 2)
            {
                scalar temp = y.c();
                y.c() = y.b();
                y.b() = y.a();
                y.a() = temp;
            }
        };

        inplaceRotate(y, triEdgei - 1 + 3);
        Swap(y.b(), y.c());
        inplaceRotate(y, 1 - triEdgej + 3);

        trii = trij;
        triEdgei = triEdgej;

        return true;
    };

    // Assume successful
    bool success = true;

    // Conform each identified edge in turn
    forAll(conformEdgeis, conformi)
    {
        const label edgei = conformEdgeis[conformi];

        // Skip if earlier operations have inadvertently resulted in this
        // edge being conformed to
        if (isConformed(edgei)) continue;

        // Get the edge and the points
        const edge e = edgePoints(edgei);
        const label pointi0 = e[0], pointi1 = e[1];

        // Find triangles which contains the first and last points
        const label patchFaceTrii0 = nextConnectedPatchFaceTri(pointi0);
        const label patchFaceTrii1 = nextConnectedPatchFaceTri(pointi1);

        // Conformation to this edge is not possible if either end point is not
        // present in the patch face triangulation
        if (patchFaceTrii0 == -1 || patchFaceTrii1 == -1) continue;

        // Track from point zero until in a triangle which contains point one.
        // Build a route of tris and edges along the track.
        DynamicList<label> routeTriis;
        DynamicList<label> routeEdgeis;
        label routePatchFaceTrii0 = -1;
        while (true)
        {
            routeTriis.clear();
            routeEdgeis.clear();

            routePatchFaceTrii0 =
                nextConnectedPatchFaceTri(pointi0, routePatchFaceTrii0 + 1);

            if (routePatchFaceTrii0 == -1) break;

            label trii = patchFaceTris[routePatchFaceTrii0];
            label triEdgei = -1;
            barycentric2D y(0, 0, 0);
            y[findIndex(triPoints(trii), pointi0)] = 1;

            forAll(patchFaceTris, iter)
            {
                if (findIndex(triPoints(trii), pointi0) != -1)
                {
                    routeTriis.clear();
                    routeEdgeis.clear();
                }

                routeTriis.append(trii);

                if (findIndex(triPoints(trii), pointi1) != -1) break;

                trackToEdge(trii, triEdgei, y, pointi1);

                if (triEdgei == -1) break;

                routeEdgeis.append(triEdges_[trii][triEdgei]);

                if (!crossEdge(trii, triEdgei, y)) break;
            }

            if (findIndex(triPoints(trii), pointi1) != -1) break;
        }

        // Conform the triangulation to the route
        if (routeEdgeis.size() == 0)
        {
            // No route was found. Don't do anything. This edge will not be
            // conformed to and the intersection will disconnect here.
            if (this->debug > 1)
            {
                WarningInFunction
                    << indent << "Failed to route edge " << e
                    << " through the triangulation of "
                    << (isSrc ? "source" : "target") << " face #"
                    << patchFacei << endl;
                writePatchFace(patchFacei, isSrc);
            }
            success = false;
        }
        else if (routeEdgeis.size() == 1)
        {
            // The route only spans a single edge. Flip it to conform.
            flipEdge(routeEdgeis.first());
        }
        else
        {
            // The route spans multiple edges. Remove all triangles in the
            // route and then split the resulting polygon along the conform
            // edge and triangulate both sides independently.

            // !!! Check that removal of the tris does not result in any points
            // becoming disconnected from the mesh
            // ...

            // Form the polygons on either side of the conform edge
            DynamicList<label> leftPolyPointis, rightPolyPointis;
            DynamicList<label> leftPolyEdgeis, rightPolyEdgeis;
            DynamicList<label> leftPolyTriis, rightPolyTriis;
            {
                // First point
                leftPolyPointis.append(pointi0);
                rightPolyPointis.append(pointi0);

                // First triangle
                {
                    const label trii = routeTriis.first();
                    const label edgei1 = routeEdgeis.first();
                    const label triEdgei1 = findIndex(triEdges_[trii], edgei1);
                    const label triEdgeiLeft = (triEdgei1 + 2) % 3;
                    const label triEdgeiRight = (triEdgei1 + 1) % 3;
                    leftPolyEdgeis.append(triEdges_[trii][triEdgeiLeft]);
                    rightPolyEdgeis.append(triEdges_[trii][triEdgeiRight]);
                    leftPolyTriis.append(trii);
                    rightPolyTriis.append(trii);
                    leftPolyPointis.append(triPoint(trii, triEdgei1));
                    rightPolyPointis.append(triPoint(trii, triEdgeiRight));
                }

                // Intermediate triangles
                for
                (
                    label routeTrii = 1;
                    routeTrii < routeTriis.size() - 1;
                    ++ routeTrii
                )
                {
                    const label trii = routeTriis[routeTrii];
                    const label edgei0 = routeEdgeis[routeTrii - 1];
                    const label edgei1 = routeEdgeis[routeTrii];
                    const label triEdgei0 = findIndex(triEdges_[trii], edgei0);
                    const label triEdgei1 = findIndex(triEdges_[trii], edgei1);
                    if ((triEdgei0 + 2) % 3 == triEdgei1)
                    {
                        const label triEdgeiLeft = (triEdgei1 + 2) % 3;
                        leftPolyEdgeis.append(triEdges_[trii][triEdgeiLeft]);
                        leftPolyTriis.append(trii);
                        leftPolyPointis.append(triPoint(trii, triEdgei1));
                    }
                    else
                    {
                        const label triEdgeiRight = (triEdgei1 + 1) % 3;
                        rightPolyEdgeis.append(triEdges_[trii][triEdgeiRight]);
                        rightPolyTriis.append(trii);
                        rightPolyPointis.append(triPoint(trii, triEdgeiRight));
                    }
                }

                // Last triangle
                {
                    const label trii = routeTriis.last();
                    const label edgei0 = routeEdgeis.last();
                    const label triEdgei0 = findIndex(triEdges_[trii], edgei0);
                    const label triEdgeiLeft = (triEdgei0 + 1) % 3;
                    const label triEdgeiRight = (triEdgei0 + 2) % 3;
                    leftPolyEdgeis.append(triEdges_[trii][triEdgeiLeft]);
                    rightPolyEdgeis.append(triEdges_[trii][triEdgeiRight]);
                    leftPolyTriis.append(trii);
                    rightPolyTriis.append(trii);
                    leftPolyPointis.append(pointi1);
                    rightPolyPointis.append(pointi1);
                }
            }

            // Disconnect and remove existing tris
            auto disconnect = [&]
            (
                const labelUList& edgeis,
                const labelList& triis
            )
            {
                forAll(edgeis, i)
                {
                    const label edgei = edgeis[i];
                    const label trii = triis[i];

                    const label edgeTrii = edgeTris_[edgei][0] != trii;
                    const label triEdgei = findIndex(triEdges_[trii], edgei);

                    const label tempEdgei = newEdgei();
                    triEdges_[trii][triEdgei] = tempEdgei;
                    edgeTris_[tempEdgei][0] = trii;
                    edgeTris_[edgei][edgeTrii] = -1;
                }
            };
            disconnect(leftPolyEdgeis, leftPolyTriis);
            disconnect(rightPolyEdgeis, rightPolyTriis);
            forAll(routeTriis, routeTrii)
            {
                removeTri(routeTriis[routeTrii]);
            }

            // Create the new conformed edge
            const label middleEdgei = newEdgei();
            leftPolyEdgeis.append(middleEdgei);
            rightPolyEdgeis.append(middleEdgei);

            // Reverse the right polygon
            SubList<label> subRightPolyPointis
            (
                rightPolyPointis,
                rightPolyPointis.size() - 1,
                1
            );
            inplaceReverseList(subRightPolyPointis);
            inplaceReverseList(rightPolyEdgeis);
            inplaceReverseList(rightPolyTriis);

            // Insert the polygons into the triangulation
            auto insert = [&]
            (
                DynamicList<label>& polyPointis,
                DynamicList<label>& polyEdgeis
            )
            {
                // Generate new edges
                polyEdgeis.resize(2*polyPointis.size() - 3, -1);
                for (label i = polyPointis.size(); i < polyEdgeis.size(); ++ i)
                {
                    polyEdgeis[i] = newEdgei();
                }

                // Triangulate
                polygonTriangulate_.triangulate
                (
                    UIndirectList<point>(points, polyPointis)
                );

                // Insert the triangulation
                forAll(polygonTriangulate_.triPoints(), trii)
                {
                    addTri
                    (
                        polygonTriangulate_.triPoints(trii, polyPointis),
                        polygonTriangulate_.triEdges(trii, polyEdgeis),
                        patchFacei,
                        isSrc
                    );
                }
            };
            insert(leftPolyPointis, leftPolyEdgeis);
            insert(rightPolyPointis, rightPolyEdgeis);
        }
    }

    // Check
    checkPatchFace(patchFacei, isSrc);

    return success;
}


template<class SrcPatchType, class TgtPatchType>
bool
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::conformPatchFaceTris
(
    const label srcFacei,
    const label tgtFacei
)
{
    return
        conformPatchFaceTris(srcFacei, tgtFacei, true)
     && conformPatchFaceTris(tgtFacei, srcFacei, false);
}


template<class SrcPatchType, class TgtPatchType>
bool
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::combinePatchFaceTris
(
    const label srcFacei,
    const label tgtFacei
)
{
    // Function to determine whether or not a given triangle is fully
    // intersected with the opposite side
    auto triIsIntersected = [&]
    (
        const label trii,
        const label otherPatchFacei,
        const triFace& otherPatchFacePatchPoints,
        const triFace& otherPatchFacePatchEdges,
        const bool isSrc
    )
    {
        const labelList& pointOtherPatchPoints =
            isSrc ? this->pointTgtPoints_ : this->pointSrcPoints_;
        const labelList& pointOtherPatchEdges =
            isSrc ? this->pointTgtEdges_ : this->pointSrcEdges_;
        const labelList& pointOtherPatchFaces =
            isSrc ? this->pointTgtFaces_ : this->pointSrcFaces_;

        forAll(triPoints_[trii], triPointi)
        {
            const label pointi = triPoint(trii, triPointi);

            if
            (
                // If the point does not intersect an other patch point ...
                findIndex
                (
                    otherPatchFacePatchPoints,
                    pointOtherPatchPoints[pointi]
                ) == -1
                // ... and does not intersect an other patch edge ...
             && findIndex
                (
                    otherPatchFacePatchEdges,
                    pointOtherPatchEdges[pointi]
                ) == -1
                // ... and not intersect the other patch face ...
             && pointOtherPatchFaces[pointi] != otherPatchFacei
            )
            {
                // ... then this triangle does not entirely intersect the
                // other patch face.
                return false;
            }
        }

        return true;
    };

    // Get the polygons which comprise all of the intersected triangles
    auto getIntersectionPolygon = [&]
    (
        const label patchFacei,
        const label otherPatchFacei,
        const bool isSrc,
        DynamicList<label>& triis,
        DynamicList<label>& edgeis
    )
    {
        const labelList& patchFaceTris =
            isSrc ? srcFaceTris_[patchFacei] : tgtFaceTris_[patchFacei];
        const labelList& triPatchFaces =
            isSrc ? triSrcFace_ : triTgtFace_;

        const triFace otherPatchFacePatchPoints
        (
            patchFacePatchPoints(otherPatchFacei, !isSrc)
        );
        const triFace otherPatchFacePatchEdges
        (
            patchFacePatchEdges(otherPatchFacei, !isSrc)
        );

        // Get an initial intersected tri
        label trii0 = -1;
        forAll(patchFaceTris, patchFaceTrii)
        {
            const label trii = patchFaceTris[patchFaceTrii];
            if
            (
                triIsIntersected
                (
                    trii,
                    otherPatchFacei,
                    otherPatchFacePatchPoints,
                    otherPatchFacePatchEdges,
                    isSrc
                )
            )
            {
                trii0 = trii;
                break;
            }
        }
        if (trii0 == -1) return;

        // Populate the star
        star::context starContext = star_.populate
        (
            trii0,
            true,
            [&](const label edgei, const label trii)
            {
                return
                    triPatchFaces[trii] == patchFacei
                 && triIsIntersected
                    (
                        trii,
                        otherPatchFacei,
                        otherPatchFacePatchPoints,
                        otherPatchFacePatchEdges,
                        isSrc
                    );
            },
            triEdges_,
            edgeTris_
        );

        // Add triangles to the list
        forAllStarFaces(star_, starTrii, trii)
        {
            triis.append(trii);
        }

        // Walk around the star polygon, adding edges to the list
        forAllStarEdges(star_, i, starEdgei, edgei)
        {
            edgeis.append(edgei);
        }
    };
    DynamicList<label> srcPolyTris, srcPolyEdges;
    DynamicList<label> tgtPolyTris, tgtPolyEdges;
    getIntersectionPolygon
    (
        srcFacei,
        tgtFacei,
        true,
        srcPolyTris,
        srcPolyEdges
    );
    getIntersectionPolygon
    (
        tgtFacei,
        srcFacei,
        false,
        tgtPolyTris,
        tgtPolyEdges
    );

    // Return if there is nothing to do
    if (srcPolyEdges.size() == 0 && tgtPolyEdges.size() == 0)
    {
        return true;
    }

    // Align the source and target polygon edges
    if (srcPolyEdges.size() && tgtPolyEdges.size())
    {
        inplaceReverseList(tgtPolyEdges);

        label tgtPolyEdgei0 = 0;
        forAll(tgtPolyEdges, tgtPolyEdgei)
        {
            const edge srcE = edgePoints(srcPolyEdges[0]);
            const edge tgtE = edgePoints(tgtPolyEdges[tgtPolyEdgei]);

            if (edge::compare(srcE, tgtE) != 0)
            {
                tgtPolyEdgei0 = tgtPolyEdgei;
                break;
            }
        }

        inplaceRotateList
        (
            static_cast<labelList&>(tgtPolyEdges),
            - tgtPolyEdgei0
        );
    }
    bool aligned;
    if (srcPolyEdges.size() == tgtPolyEdges.size())
    {
        aligned = true;
        forAll(srcPolyEdges, polyEdgei)
        {
            const edge srcE = edgePoints(srcPolyEdges[polyEdgei]);
            const edge tgtE = edgePoints(tgtPolyEdges[polyEdgei]);

            if (edge::compare(srcE, tgtE) == 0)
            {
                aligned = false;
                break;
            }
        }
    }
    else
    {
        aligned = false;
    }
    if (!aligned)
    {
        if (this->debug > 1)
        {
            WarningInFunction
                << indent << "Failed to combine intersected parts of source "
                << "face #" << srcFacei << " and target face #" << tgtFacei
                << endl;
            writePatchFace(srcFacei, true);
            writePatchFace(tgtFacei, false);
        }
        return false;
    }

    // Mark triangles that belong to either polygon as candidates
    forAll(srcPolyTris, srcPolyTrii)
    {
        const label srcTrii = srcPolyTris[srcPolyTrii];
        markedTriTris_.append(srcTrii);
        triMarkedTris_[srcTrii] = markedTriTris_.size() - 1;
    }
    forAll(tgtPolyTris, tgtPolyTrii)
    {
        const label tgtTrii = tgtPolyTris[tgtPolyTrii];
        markedTriTris_.append(tgtTrii);
        triMarkedTris_[tgtTrii] = markedTriTris_.size() - 1;
    }

    // Change the aligned edges so that they connect opposite sides, thereby
    // disconnecting the intersected triangles from the rest of their surfaces.
    // Tris in the polygon are connected to the source edge and non-star tris
    // are connected to the target edge. If there are no non-star tris
    // connected to a pair of aligned edges then the target edge will not be
    // connected to any triangles as a result of this operation. However, we do
    // not remove the target edge in this case because we need it later to
    // attach the new intersected face to.
    forAll(srcPolyEdges, polyEdgei)
    {
        const label srcEdgei = srcPolyEdges[polyEdgei];
        const label tgtEdgei = tgtPolyEdges[polyEdgei];

        if (srcEdgei == tgtEdgei) continue;

        const label srcEdgeTrii =
            edgeTris_[srcEdgei][0] == -1
         || triMarkedTris_[edgeTris_[srcEdgei][0]] == -1
         || triSrcFace_[edgeTris_[srcEdgei][0]] == -1;
        const label tgtEdgeTrii =
            edgeTris_[tgtEdgei][0] == -1
         || triMarkedTris_[edgeTris_[tgtEdgei][0]] == -1
         || triTgtFace_[edgeTris_[tgtEdgei][0]] == -1;

        const label tgtTrii = edgeTris_[tgtEdgei][tgtEdgeTrii];
        const label srcTrij = edgeTris_[srcEdgei][!srcEdgeTrii];
        const label tgtTrij = edgeTris_[tgtEdgei][!tgtEdgeTrii];

        edgeTris_[srcEdgei][!srcEdgeTrii] = tgtTrii;
        edgeTris_[tgtEdgei][tgtEdgeTrii] = srcTrij;

        const label tgtTriEdgei = findIndex(triEdges_[tgtTrii], tgtEdgei);
        triEdges_[tgtTrii][tgtTriEdgei] = srcEdgei;

        if (srcTrij != -1)
        {
            const label srcTriEdgej = findIndex(triEdges_[srcTrij], srcEdgei);
            triEdges_[srcTrij][srcTriEdgej] = tgtEdgei;
        }

        if (srcTrij != -1 && tgtTrij != -1)
        {
            edgeFrontEdges_[tgtEdgei] = frontEdgeEdges_.size();
            frontEdgeEdges_.append(tgtEdgei);
        }
    }

    // Circle the edges to create the intersection face
    const label facei = this->faces().size();
    this->faces_.append(face(srcPolyEdges.size()));
    faceEdges_.append(labelList(srcPolyEdges.size()));
    this->srcFaceFaces_[srcFacei].append(facei);
    this->tgtFaceFaces_[tgtFacei].append(facei);
    this->faceSrcFaces_.append(srcFacei);
    this->faceTgtFaces_.append(tgtFacei);
    forAll(srcPolyEdges, polyEdgei)
    {
        const label srcEdgei = srcPolyEdges[polyEdgei];
        const label srcEdgeTrii =
            edgeTris_[srcEdgei][0] == -1
         || triMarkedTris_[edgeTris_[srcEdgei][0]] == -1
         || triSrcFace_[edgeTris_[srcEdgei][0]] == -1;

        const label srcTrii = edgeTris_[srcEdgei][srcEdgeTrii];
        const label srcTriEdgei = findIndex(triEdges_[srcTrii], srcEdgei);

        this->faces_.last()[polyEdgei] =
            triEdgePoints(srcTrii, srcTriEdgei).start();

        const label tgtEdgei = tgtPolyEdges[polyEdgei];

        faceEdges_.last()[polyEdgei] = tgtEdgei;
        intersectEdgeFaces_[tgtEdgei][intersectEdgeFaces_[tgtEdgei][0] != -1] =
            facei;
    }

    // Clear the marked tris
    forAll(markedTriTris_, candidateTrii)
    {
        const label trii = markedTriTris_[candidateTrii];
        if (trii != -1)
        {
            triMarkedTris_[trii] = -1;
        }
    }
    markedTriTris_.clear();

    // Check
    checkPatchFace(srcFacei, true);
    checkPatchFace(tgtFacei, false);

    // Remove the intersection triangles
    forAll(srcPolyTris, srcPolyTrii)
    {
        removeTri(srcPolyTris[srcPolyTrii]);
    }
    forAll(tgtPolyTris, tgtPolyTrii)
    {
        removeTri(tgtPolyTris[tgtPolyTrii]);
    }

    // Check
    checkPatchFace(srcFacei, true);
    checkPatchFace(tgtFacei, false);

    // Update the propagation front
    forAll(tgtPolyEdges, polyEdgei)
    {
        const label edgei = tgtPolyEdges[polyEdgei];

        const label trii0 = edgeTris_[edgei][0];
        const label trii1 = edgeTris_[edgei][1];

        const bool isFront =
            trii0 != -1
         && trii1 != -1
         && (triSrcFace_[trii0] == -1) != (triSrcFace_[trii1] == -1);

        if (isFront && edgeFrontEdges_[edgei] == -1)
        {
            edgeFrontEdges_[edgei] = frontEdgeEdges_.size();
            frontEdgeEdges_.append(edgei);
        }

        if (!isFront && edgeFrontEdges_[edgei] != -1)
        {
            frontEdgeEdges_[edgeFrontEdges_[edgei]] = -1;
            edgeFrontEdges_[edgei] = -1;
        }
    }

    // Check
    checkPatchFace(srcFacei, true);
    checkPatchFace(tgtFacei, false);

    return true;
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::initialise
(
    const vectorField& srcPointNormals
)
{
    // Clear the base class data ...

    this->points_.clear();

    this->srcPointPoints_ = -1;
    this->tgtPointPoints_ = -1;
    this->pointSrcPoints_.clear();
    this->pointTgtPoints_.clear();

    forAll(this->srcEdgePoints_, srcEdgei)
    {
        this->srcEdgePoints_[srcEdgei].clear();
    }
    forAll(this->tgtEdgePoints_, tgtEdgei)
    {
        this->tgtEdgePoints_[tgtEdgei].clear();
    }
    this->pointSrcEdges_.clear();
    this->pointTgtEdges_.clear();

    this->pointSrcFaces_.clear();
    this->pointTgtFaces_.clear();

    this->faces_.clear();

    forAll(this->srcFaceFaces_, srcFacei)
    {
        this->srcFaceFaces_[srcFacei].clear();
    }
    forAll(this->tgtFaceFaces_, tgtFacei)
    {
        this->tgtFaceFaces_[tgtFacei].clear();
    }
    this->faceSrcFaces_.clear();
    this->faceTgtFaces_.clear();

    // Initialise with the source and target patch data ...

    // Points
    const label nPoints =
        this->srcPatch_.nPoints() + this->tgtPatch_.nPoints();
    srcPoints_.resize(nPoints, point::uniform(NaN));
    srcPointNormals_.resize(nPoints, vector::uniform(NaN));
    tgtPoints_.resize(nPoints, point::uniform(NaN));
    pointPoints_.resize(nPoints, -1);
    this->pointSrcPoints_.resize(nPoints, -1);
    this->pointTgtPoints_.resize(nPoints, -1);
    this->pointSrcEdges_.resize(nPoints, -1);
    this->pointTgtEdges_.resize(nPoints, -1);
    this->pointSrcFaces_.resize(nPoints, -1);
    this->pointTgtFaces_.resize(nPoints, -1);
    forAll(this->srcPatch_.localPoints(), srcPointi)
    {
        const label pointi = srcPointi;

        srcPoints_[pointi] = this->srcPatch_.localPoints()[srcPointi];
        tgtPoints_[pointi] = this->srcPatch_.localPoints()[srcPointi];
        srcPointNormals_[pointi] = srcPointNormals[srcPointi];
        pointPoints_[pointi] = pointi;
        this->srcPointPoints_[srcPointi] = pointi;
        this->pointSrcPoints_[pointi] = srcPointi;
    }
    forAll(this->tgtPatch_.localPoints(), tgtPointi)
    {
        const label pointi = this->srcPatch_.nPoints() + tgtPointi;

        srcPoints_[pointi] = this->tgtPatch_.localPoints()[tgtPointi];
        tgtPoints_[pointi] = this->tgtPatch_.localPoints()[tgtPointi];
        pointPoints_[pointi] = pointi;
        this->tgtPointPoints_[tgtPointi] = pointi;
        this->pointTgtPoints_[pointi] = tgtPointi;
    }

    // Edges
    const label nEdges = this->srcPatch_.nEdges() + this->tgtPatch_.nEdges();
    edgeTris_.resize(nEdges, labelPair(-1, -1));
    intersectEdgeFaces_.resize(nEdges, labelPair(-1, -1));
    forAll(this->srcPatch_.faceEdges(), srcFacei)
    {
        const label trii = srcFacei;

        forAll(this->srcPatch_.faceEdges()[srcFacei], srcFaceEdgei)
        {
            const label srcEdgei =
                this->srcPatch_.faceEdges()[srcFacei][srcFaceEdgei];
            const label edgei = srcEdgei;

            const edge& e = this->srcPatch_.edges()[srcEdgei];
            const edge fe =
                this->srcPatch_.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

            edgeTris_[edgei][edge::compare(e, fe) < 0] = trii;
        }
    }
    forAll(this->srcPatch_.edges(), srcEdgei)
    {
        const edge& e = this->srcPatch_.edges()[srcEdgei];

        this->srcEdgePoints_[srcEdgei].append(this->srcPointPoints_[e[0]]);
        this->srcEdgePoints_[srcEdgei].append(this->srcPointPoints_[e[1]]);
    }
    forAll(this->tgtPatch_.faceEdges(), tgtFacei)
    {
        const label trii = this->srcPatch_.size() + tgtFacei;

        forAll(this->tgtPatch_.faceEdges()[tgtFacei], tgtFaceEdgei)
        {
            const label tgtEdgei =
                this->tgtPatch_.faceEdges()[tgtFacei][tgtFaceEdgei];
            const label edgei = this->srcPatch_.nEdges() + tgtEdgei;

            const edge& e = this->tgtPatch_.edges()[tgtEdgei];
            const edge fe =
                this->tgtPatch_.localFaces()[tgtFacei].faceEdge(tgtFaceEdgei);

            edgeTris_[edgei][edge::compare(e, fe) < 0] = trii;
        }
    }
    forAll(this->tgtPatch_.edges(), tgtEdgei)
    {
        const edge& e = this->tgtPatch_.edges()[tgtEdgei];

        this->tgtEdgePoints_[tgtEdgei].append(this->tgtPointPoints_[e[0]]);
        this->tgtEdgePoints_[tgtEdgei].append(this->tgtPointPoints_[e[1]]);
    }

    // Tris
    const label nTris = this->srcPatch_.size() + this->tgtPatch_.size();
    triPoints_.resize(nTris);
    triEdges_.resize(nTris);
    triSrcFace_.resize(nTris, -1);
    triTgtFace_.resize(nTris, -1);
    forAll(this->srcPatch_.localFaces(), srcFacei)
    {
        const label trii = srcFacei;

        forAll(this->srcPatch_.localFaces()[srcFacei], i)
        {
            triPoints_[trii][i] = this->srcPatch_.localFaces()[srcFacei][i];
            triEdges_[trii][i] = this->srcPatch_.faceEdges()[srcFacei][i];
        }
        triSrcFace_[trii] = srcFacei;
        srcFaceTris_[srcFacei].resize(1, trii);
    }
    forAll(this->tgtPatch_.localFaces(), tgtFacei)
    {
        const label trii = this->srcPatch_.size() + tgtFacei;

        forAll(this->tgtPatch_.localFaces()[tgtFacei], i)
        {
            triPoints_[trii][i] =
                this->srcPatch_.nPoints()
              + this->tgtPatch_.localFaces()[tgtFacei][i];
            triEdges_[trii][i] =
                this->srcPatch_.nEdges()
              + this->tgtPatch_.faceEdges()[tgtFacei][i];
        }
        triTgtFace_[trii] = tgtFacei;
        tgtFaceTris_[tgtFacei].resize(1, trii);
    }

    // Removal
    removedEdges_.clear();
    removedTris_.clear();

    // Front propagation
    frontEdgeEdges_.clear();
    edgeFrontEdges_ = DynamicList<label>(edgeTris_.size(), -1);

    // Insertion
    candidateTriTris_.clear();
    triCandidateTris_ = DynamicList<label>(triPoints_.size(), -1);

    // Marking
    markedTriTris_.clear();
    triMarkedTris_ = DynamicList<label>(triPoints_.size(), -1);

    checkPatchFaces(true);
    checkPatchFaces(false);
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::clean()
{
    // Resolve point-points
    inplaceRenumber(pointPoints_, triPoints_);
    inplaceRenumber(pointPoints_, this->srcPointPoints_);
    inplaceRenumber(pointPoints_, this->tgtPointPoints_);
    inplaceRenumber(pointPoints_, this->srcEdgePoints_);
    inplaceRenumber(pointPoints_, this->tgtEdgePoints_);

    // Remove points by shuffling up
    labelList oldPointNewPoints(pointPoints_.size(), -1);
    label pointi = 0;
    forAll(pointPoints_, pointj)
    {
        if (pointPoints_[pointj] == pointj)
        {
            oldPointNewPoints[pointj] = pointi;
            srcPoints_[pointi] = srcPoints_[pointj];
            srcPointNormals_[pointi] = srcPointNormals_[pointj];
            tgtPoints_[pointi] = tgtPoints_[pointj];
            pointPoints_[pointi] = pointi;
            this->pointSrcPoints_[pointi] = this->pointSrcPoints_[pointj];
            this->pointTgtPoints_[pointi] = this->pointTgtPoints_[pointj];
            this->pointSrcEdges_[pointi] = this->pointSrcEdges_[pointj];
            this->pointTgtEdges_[pointi] = this->pointTgtEdges_[pointj];
            this->pointSrcFaces_[pointi] = this->pointSrcFaces_[pointj];
            this->pointTgtFaces_[pointi] = this->pointTgtFaces_[pointj];
            ++ pointi;
        }
    }
    srcPoints_.resize(pointi);
    srcPointNormals_.resize(pointi);
    tgtPoints_.resize(pointi);
    pointPoints_.resize(pointi);
    this->pointSrcPoints_.resize(pointi);
    this->pointTgtPoints_.resize(pointi);
    this->pointSrcEdges_.resize(pointi);
    this->pointTgtEdges_.resize(pointi);
    this->pointSrcFaces_.resize(pointi);
    this->pointTgtFaces_.resize(pointi);

    // Remove edges by shuffling up
    labelList oldEdgeNewEdges(edgeTris_.size(), -1);
    label edgei = 0;
    forAll(edgeTris_, edgej)
    {
        if
        (
            edgeTris_[edgej] != labelPair(-1, -1)
         || intersectEdgeFaces_[edgej] != labelPair(-1, -1)
        )
        {
            oldEdgeNewEdges[edgej] = edgei;
            edgeTris_[edgei] = edgeTris_[edgej];
            intersectEdgeFaces_[edgei] = intersectEdgeFaces_[edgej];
            edgeFrontEdges_[edgei] = edgeFrontEdges_[edgej];
            ++ edgei;
        }
    }
    edgeTris_.resize(edgei);
    intersectEdgeFaces_.resize(edgei);
    edgeFrontEdges_.resize(edgei);

    // Remove tris by shuffling up
    labelList oldTriNewTris(triPoints_.size(), -1);
    label trii = 0;
    forAll(triPoints_, trij)
    {
        if (triPoints_[trij] != FixedList<label, 3>({-1, -1, -1}))
        {
            oldTriNewTris[trij] = trii;
            triPoints_[trii] = triPoints_[trij];
            triEdges_[trii] = triEdges_[trij];
            triSrcFace_[trii] = triSrcFace_[trij];
            triTgtFace_[trii] = triTgtFace_[trij];
            ++ trii;
        }
    }
    triPoints_.resize(trii);
    triEdges_.resize(trii);
    triSrcFace_.resize(trii);
    triTgtFace_.resize(trii);

    // Map
    inplaceRenumber(oldTriNewTris, edgeTris_);
    inplaceRenumber(oldPointNewPoints, triPoints_);
    inplaceRenumber(oldEdgeNewEdges, triEdges_);
    inplaceRenumber(oldTriNewTris, srcFaceTris_);
    inplaceRenumber(oldTriNewTris, tgtFaceTris_);
    inplaceRenumber(oldEdgeNewEdges, frontEdgeEdges_);
    inplaceRenumber(oldPointNewPoints, this->faces_);
    inplaceRenumber(oldEdgeNewEdges, faceEdges_);
    inplaceRenumber(oldPointNewPoints, this->srcPointPoints_);
    inplaceRenumber(oldPointNewPoints, this->tgtPointPoints_);
    inplaceRenumber(oldPointNewPoints, this->srcEdgePoints_);
    inplaceRenumber(oldPointNewPoints, this->tgtEdgePoints_);

    // Removal
    removedEdges_.clear();
    removedTris_.clear();

    // Insertion
    candidateTriTris_.clear();
    triCandidateTris_ = DynamicList<label>(triPoints_.size(), -1);

    // Marking
    markedTriTris_.clear();
    triMarkedTris_ = DynamicList<label>(triPoints_.size(), -1);

    checkPatchFaces(true);
    checkPatchFaces(false);
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::finalise()
{
    clean();

    // Add the remaining triangles as uncoupled faces ...

    const label nFaces = this->faces_.size();

    this->faces_.resize(nFaces + triPoints_.size());
    faceEdges_.resize(nFaces + triPoints_.size());
    this->faceSrcFaces_.resize(nFaces + triPoints_.size());
    this->faceTgtFaces_.resize(nFaces + triPoints_.size());

    forAll(triPoints_, trii)
    {
        const label facei = nFaces + trii;

        if (triSrcFace_[trii] != -1)
        {
            this->faces_[facei] = triPoints(trii);
            faceEdges_[facei] = labelList(triEdges_[trii]);
            this->srcFaceFaces_[triSrcFace_[trii]].append(facei);
            this->faceSrcFaces_[facei] = triSrcFace_[trii];
            this->faceTgtFaces_[facei] = -1;
        }
        else
        {
            this->faces_[facei] = triPoints(trii).reverseFace();
            faceEdges_[facei] = labelList(reverseList(triEdges_[trii]));
            this->tgtFaceFaces_[triTgtFace_[trii]].append(facei);
            this->faceSrcFaces_[facei] = -1;
            this->faceTgtFaces_[facei] = triTgtFace_[trii];
        }
    }

    nonIntersectEdgeFaces_.resize(edgeTris_.size());
    forAll(edgeTris_, edgei)
    {
        forAll(edgeTris_[edgei], edgeTrii)
        {
            nonIntersectEdgeFaces_[edgei][edgeTrii] =
                edgeTris_[edgei][edgeTrii] + nFaces;
        }
    }

    checkPatchFaces(true);
    checkPatchFaces(false);
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::unFinalise()
{
    clean();

    // Remove all uncoupled faces and store them as triangles ...

    label nFaces = 0;
    while
    (
        nFaces < this->faces_.size()
     && this->faceSrcFaces_[nFaces] != -1
     && this->faceTgtFaces_[nFaces] != -1
    )
    {
        ++ nFaces;
    }

    this->faces_.resize(nFaces);
    faceEdges_.resize(nFaces);
    this->faceSrcFaces_.resize(nFaces);
    this->faceTgtFaces_.resize(nFaces);

    forAll(triPoints_, trii)
    {
        DynamicList<label>& faceis =
            triSrcFace_[trii] != -1
          ? this->srcFaceFaces_[triSrcFace_[trii]]
          : this->tgtFaceFaces_[triTgtFace_[trii]];

        label n = 0;
        while (n < faceis.size() && faceis[n] < nFaces)
        {
            ++ n;
        }

        faceis.resize(n);
    }

    nonIntersectEdgeFaces_.clear();

    checkPatchFaces(true);
    checkPatchFaces(false);
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::write()
{
    if (this->debug > 2)
    {
        finalise();

        // Use base class to write the patch
        this->report(name(writei_));

        unFinalise();

        // Write the edge front
        const fileName frontFileName =
            type() + "_front_" + name(writei_) + ".vtk";
        Info<< indent << "Writing front to " << frontFileName << endl;

        DynamicList<point> writePoints(frontEdgeEdges_.size()*2);
        DynamicList<labelPair> writeLines(frontEdgeEdges_.size()*2);
        forAll(frontEdgeEdges_, frontEdgei)
        {
            const label edgei = frontEdgeEdges_[frontEdgei];

            if (edgei == -1) continue;

            const edge e = edgePoints(edgei);
            writePoints.append(tgtPoints_[e.start()]);
            writePoints.append(tgtPoints_[e.end()]);

            const label i = writePoints.size() - 2;
            writeLines.append(labelPair(i, i + 1));
        }

        vtkWritePolyData::write
        (
            frontFileName,
            "front",
            false,
            writePoints,
            labelList(),
            writeLines,
            faceList()
        );

        writei_ ++;
    }
}


template<class SrcPatchType, class TgtPatchType>
void Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::writePatchFace
(
    const label patchFacei,
    const bool isSrc
) const
{
    OFstream os
    (
        word(isSrc ? "src" : "tgt") + "Face_" + name(patchFacei) + ".obj"
    );

    OFstream tos
    (
        word(isSrc ? "src" : "tgt") + "FaceTris_" + name(patchFacei) + ".obj"
    );

    Info<< indent << "Writing patch face to " << os.name()
        << " and patch face triangulation to " << tos.name()
        << incrIndent << endl;

    forAll(patchFacePatchPoints(patchFacei, isSrc), patchFacePatchPointi)
    {
        const label patchPointi =
            patchFacePatchPoints(patchFacei, isSrc)[patchFacePatchPointi];
        const point& p =
            isSrc
          ? this->srcPatch_.localPoints()[patchPointi]
          : this->tgtPatch_.localPoints()[patchPointi];
        os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }
    os  << "f 1 2 3" << nl;

    const labelList& patchFaceTris =
        isSrc ? srcFaceTris_[patchFacei] : tgtFaceTris_[patchFacei];

    forAll(patchFaceTris, patchFaceTrii)
    {
        const label trii = patchFaceTris[patchFaceTrii];

        Info<< indent << "tri #" << trii << " points=" << triPoints(trii)
            << " edges=" << triEdges_[trii] << endl;

        forAll(triPoints_[trii], triPointi)
        {
            const label pointi = triPoint(trii, triPointi);
            const point& p = isSrc ? srcPoints_[pointi] : tgtPoints_[pointi];
            tos << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }
        tos << "f";
        forAll(triPoints_[trii], triPointi)
        {
            tos << " " << 1 + 3*patchFaceTrii + triPointi;
        }
        tos << nl;
    }

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SrcPatchType, class TgtPatchType>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::TriPatchIntersection
(
    const SrcPatchType& srcPatch,
    const TgtPatchType& tgtPatch,
    const scalar snapTol
)
:
    TriPatchIntersection(srcPatch, srcPatch.pointNormals(), tgtPatch, snapTol)
{}


template<class SrcPatchType, class TgtPatchType>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::TriPatchIntersection
(
    const SrcPatchType& srcPatch,
    const vectorField& srcPointNormals,
    const TgtPatchType& tgtPatch,
    const scalar snapTol
)
:
    PatchIntersection<SrcPatchType, TgtPatchType>(srcPatch, tgtPatch),

    srcPoints_(),
    srcPointNormals_(),
    tgtPoints_(this->points_),
    pointPoints_(),

    edgeTris_(),
    intersectEdgeFaces_(),
    nonIntersectEdgeFaces_(),

    triPoints_(),
    triEdges_(),
    triSrcFace_(),
    triTgtFace_(),
    srcFaceTris_(srcPatch.size()),
    tgtFaceTris_(tgtPatch.size()),

    faceEdges_(),

    removedTris_(),
    removedEdges_(),

    frontEdgeEdges_(),
    edgeFrontEdges_(),

    candidateTriTris_(),
    triCandidateTris_(),

    markedTriTris_(),
    triMarkedTris_(),

    polygonTriangulate_(),

    star_(),

    writei_(0)
{
    cpuTime time;

    Info<< indent << type() << ": Intersecting "
        << this->srcPatch_.size() << " source tri faces and "
        << this->tgtPatch_.size() << " target tri faces" << incrIndent << endl;

    if (this->debug)
    {
        Info<< indent << "Writing tri patches" << incrIndent << endl;
        const fileName srcFileName = type() + "_srcPatch.vtk";
        Info<< indent << "Writing patch to " << srcFileName << endl;
        vtkWritePolyData::write
        (
            srcFileName,
            "source",
            false,
            this->srcPatch_.localPoints(),
            labelList(),
            labelListList(),
            this->srcPatch_.localFaces(),
            "normals",
            true,
            srcPointNormals
        );
        const fileName tgtFileName = type() + "_tgtPatch.vtk";
        Info<< indent << "Writing patch to " << tgtFileName << endl;
        vtkWritePolyData::write
        (
            tgtFileName,
            "target",
            false,
            this->tgtPatch_.localPoints(),
            labelList(),
            labelListList(),
            this->tgtPatch_.localFaces()
        );
        Info<< decrIndent;
    }

    // Create bound spheres for patch faces for proximity testing
    List<Tuple2<point, scalar>> srcFaceSpheres(this->srcPatch_.size());
    List<Tuple2<point, scalar>> tgtFaceSpheres(this->tgtPatch_.size());
    forAll(this->srcPatch_, srcFacei)
    {
        const triFace& srcFace = this->srcPatch_.localFaces()[srcFacei];
        srcFaceSpheres[srcFacei] =
            trivialBoundSphere
            (
                this->srcPatch_.localPoints(),
                {srcFace[0], srcFace[1], srcFace[2], -1},
                3
            );
        srcFaceSpheres[srcFacei].second() *= 1 + snapTol;
    }
    forAll(this->tgtPatch_, tgtFacei)
    {
        const triFace& tgtFace = this->tgtPatch_.localFaces()[tgtFacei];
        tgtFaceSpheres[tgtFacei] =
            trivialBoundSphere
            (
                this->tgtPatch_.localPoints(),
                {tgtFace[0], tgtFace[1], tgtFace[2], -1},
                3
            );
        tgtFaceSpheres[tgtFacei].second() *= 1 + snapTol;
    }

    // Construct table to store what faces have been snapped
    HashSet<labelPair, labelPair::Hash<>> srcFaceTgtFaceSnaps
    (
        12*(this->srcPatch_.size() + this->tgtPatch_.size())
    );

    // Count the number of successes and failures
    label nIntersections = 0, nIntersectionFailures = 0;

    // Function for intersecting two patch faces
    auto intersect = [&](const label srcFacei, const label tgtFacei)
    {
        // Single snapping stage
        /*
        snapPatchFaceTris(srcFacei, tgtFacei, snapTol);
        */

        // Propagate out and snap everything in advance of the intersections
        const List<triFace>& srcLocalFaces = this->srcPatch_.localFaces();
        const List<triFace>& tgtLocalFaces = this->tgtPatch_.localFaces();
        const labelListList& srcPointFaces = this->srcPatch_.pointFaces();
        const labelListList& tgtPointFaces = this->tgtPatch_.pointFaces();
        forAll(srcLocalFaces[srcFacei], srcFacePointi)
        {
            const label srcPointi =
                srcLocalFaces[srcFacei][srcFacePointi];

            forAll(srcPointFaces[srcPointi], srcPointFacei)
            {
                const label srcFacej =
                    srcPointFaces[srcPointi][srcPointFacei];

                forAll(tgtLocalFaces[tgtFacei], tgtFacePointi)
                {
                    const label tgtPointi =
                        tgtLocalFaces[tgtFacei][tgtFacePointi];

                    forAll(tgtPointFaces[tgtPointi], tgtPointFacei)
                    {
                        const label tgtFacej =
                            tgtPointFaces[tgtPointi][tgtPointFacei];

                        const point& srcC = srcFaceSpheres[srcFacej].first();
                        const scalar srcR = srcFaceSpheres[srcFacej].second();
                        const point& tgtC = tgtFaceSpheres[tgtFacej].first();
                        const scalar tgtR = tgtFaceSpheres[tgtFacej].second();

                        if
                        (
                            magSqr(srcC - tgtC) < sqr(srcR + tgtR)
                         && !srcFaceTgtFaceSnaps.found({srcFacej, tgtFacej})
                        )
                        {
                            srcFaceTgtFaceSnaps.insert({srcFacej, tgtFacej});
                            snapPatchFaceTris(srcFacej, tgtFacej, snapTol);
                        }
                    }
                }
            }
        }

        intersectPatchFaceTris(srcFacei, tgtFacei);

        write();

        const bool conformFailure = !conformPatchFaceTris(srcFacei, tgtFacei);

        const bool combineFailure = !combinePatchFaceTris(srcFacei, tgtFacei);

        nIntersections ++;
        nIntersectionFailures += conformFailure || combineFailure;

        write();
    };

    // Target search tree. Used to find an initial source and target face to
    // intersect. Once the first intersection has been done the rest follow via
    // a front-propagation algorithm, without reference to this tree.
    indexedOctree<treeDataPrimitivePatch<TgtPatchType>> tgtTree
    (
        treeDataPrimitivePatch<TgtPatchType>
        (
            false,
            this->tgtPatch_,
            indexedOctree<treeDataPrimitivePatch<TgtPatchType>>::perturbTol()
        ),
        treeBoundBox(this->tgtPatch_.localPoints()).extend(1e-4),
        8,
        10,
        3
    );

    // Find intersection operation, excluding a hash set of previously obtained
    // intersections. Used to get multiple target triangles intersected by a
    // single source point and point normal (i.e., two target triangles if the
    // source point intersects a target edge).
    class findIntersectExcludingOp
    {
    public:

        const indexedOctree<treeDataPrimitivePatch<TgtPatchType>>& tree_;

        const labelHashSet& excludingIndices_;

    public:

        //- Construct from components
        findIntersectExcludingOp
        (
            const indexedOctree<treeDataPrimitivePatch<TgtPatchType>>& tree,
            const labelHashSet& hitIndices
        )
        :
            tree_(tree),
            excludingIndices_(hitIndices)
        {}

        //- Calculate intersection of triangle with ray. Sets result
        //  accordingly
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const
        {
            if (excludingIndices_.found(index))
            {
                return false;
            }
            else
            {
                typename
                    treeDataPrimitivePatch<TgtPatchType>::findIntersectOp
                    iop(tree_);

                return
                    iop
                    (
                        index,
                        start,
                        end,
                        intersectionPoint
                    );
            }
        }
    };

    // Populate local data from the source and target patches
    initialise(srcPointNormals);
    write();

    // Loop the source points, looking for ones that have not been intersected
    forAll(this->points_, pointi)
    {
        // Get the next source patch point
        const label srcPointi = this->pointSrcPoints_[pointi];

        // Continue if this point is already intersected
        if
        (
            srcPointi == -1
         || this->pointTgtPoints_[pointi] != -1
         || this->pointTgtEdges_[pointi] != -1
         || this->pointTgtFaces_[pointi] != -1
        ) continue;

        // Get the projection geometry for this source point
        const point& srcP = srcPoints_[pointi];
        const vector& srcN = srcPointNormals_[pointi];
        scalar srcL = 0;
        forAll(this->srcPatch_.pointEdges()[srcPointi], srcPointEdgei)
        {
            const label srcEdgei =
                this->srcPatch_.pointEdges()[srcPointi][srcPointEdgei];
            srcL +=
                this->srcPatch_.edges()[srcEdgei].mag
                (
                    this->srcPatch_.localPoints()
                );
        }
        srcL /= this->srcPatch_.pointEdges()[srcPointi].size();

        // Find all the target faces that this source point projects to
        labelHashSet tgtFaceis;
        while (true)
        {
            pointIndexHit hit =
                tgtTree.findLine
                (
                    srcP - srcL*srcN,
                    srcP + srcL*srcN,
                    findIntersectExcludingOp(tgtTree, tgtFaceis)
                );

            if (!hit.hit()) break;

            tgtFaceis.insert(hit.index());
        }

        // Continue if this point does not project to the opposite patch
        if (tgtFaceis.empty()) continue;

        // Loop all the potential source/target face intersections until an
        // edge front is generated
        forAll(this->srcPatch_.pointFaces()[srcPointi], srcPointFacei)
        {
            const label srcFacei =
                this->srcPatch_.pointFaces()[srcPointi][srcPointFacei];

            forAllConstIter(labelHashSet, tgtFaceis, tgtFaceiIter)
            {
                intersect(srcFacei, tgtFaceiIter.key());

                if (frontEdgeEdges_.size()) break;
            }

            if (frontEdgeEdges_.size()) break;
        }

        // Propagate until the edge front is empty
        while (frontEdgeEdges_.size())
        {
            const label edgei = frontEdgeEdges_.remove();

            if (edgei == -1) continue;

            edgeFrontEdges_[edgei] = -1;

            label srcTrii = edgeTris_[edgei][0];
            label tgtTrii = edgeTris_[edgei][1];

            if (triSrcFace_[srcTrii] == -1)
            {
                Swap(srcTrii, tgtTrii);
            }

            intersect(triSrcFace_[srcTrii], triTgtFace_[tgtTrii]);
        }
    }

    // Populate data in the base class which was not generated as part of the
    // intersection process
    finalise();

    this->report();

    // Warn about any failures
    if (nIntersectionFailures)
    {
        Info<< indent << "*** Topology could not be generated in "
            << nIntersectionFailures << "/" << nIntersections << " cases"
            << endl << indent << "    The intersection may be incomplete"
            << endl;
    }

    Info<< indent << this->faces_.size() << " faces generated in "
            << time.cpuTimeIncrement() << 's' << endl;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SrcPatchType, class TgtPatchType>
Foam::TriPatchIntersection<SrcPatchType, TgtPatchType>::
~TriPatchIntersection()
{}


// ************************************************************************* //
