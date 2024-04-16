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

#include "FacePatchIntersection.H"
#include "cpuTime.H"
#include "primitiveTriPatch.H"
#include "polygonTriangulate.H"
#include "TriPatchIntersection.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T, class ListType>
void inplaceRemove(ListType& lst, const T& value)
{
    if (lst.size())
    {
        inplaceRemove(lst, value, lst[0]);
    }
}

template<class T, class NotT, class ListType>
void inplaceRemove(ListType& lst, const T& value, const NotT&)
{
    forAll(lst, i)
    {
        if (lst[i].size())
        {
            inplaceRemove(lst[i], value, lst[i][0]);
        }
    }
}

template<class T, class ListType>
void inplaceRemove(ListType& lst, const T& value, const T&)
{
    label iNew = 0;
    forAll(lst, iOld)
    {
        if (lst[iOld] != value)
        {
            lst[iNew] = lst[iOld];
            ++ iNew;
        }
    }
    lst.resize(iNew);
}

template<class ListType>
ListType reorder
(
    const label size,
    const typename ListType::value_type& defaultValue,
    const labelUList& oldToNew,
    const ListType& lst
)
{
    ListType newLst(size, defaultValue);

    forAll(lst, elemI)
    {
        if (oldToNew[elemI] >= 0)
        {
            newLst[oldToNew[elemI]] = lst[elemI];
        }
    }

    return newLst;
}

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SrcPatchType, class TgtPatchType>
Foam::FacePatchIntersection<SrcPatchType, TgtPatchType>::FacePatchIntersection
(
    const SrcPatchType& srcPatch,
    const TgtPatchType& tgtPatch,
    const scalar snapTol
)
:
    FacePatchIntersection(srcPatch, srcPatch.pointNormals(), tgtPatch, snapTol)
{}


template<class SrcPatchType, class TgtPatchType>
Foam::FacePatchIntersection<SrcPatchType, TgtPatchType>::FacePatchIntersection
(
    const SrcPatchType& srcPatch,
    const vectorField& srcPointNormals,
    const TgtPatchType& tgtPatch,
    const scalar snapTol
)
:
    PatchIntersection<SrcPatchType, TgtPatchType>(srcPatch, tgtPatch)
{
    cpuTime time;

    Info<< type() << ": Intersecting "
        << this->srcPatch_.size() << " source faces and "
        << this->tgtPatch_.size() << " target faces" << incrIndent << endl;

    if (this->debug)
    {
        Info<< indent << "Writing patches" << incrIndent << endl;
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

    // Step 1: Triangulate the source and target patches and create addressing
    DynamicList<triFace> srcTris, tgtTris;
    DynamicList<label> srcTriSrcFaces, tgtTriTgtFaces;
    DynamicList<FixedList<label, 3>> srcTriSrcEdges, tgtTriTgtEdges;
    polygonTriangulate triEngine;
    forAll(this->srcPatch_, srcFacei)
    {
        const face& srcFacePoints = this->srcPatch_[srcFacei];
        const labelList& srcFaceEdges = this->srcPatch_.faceEdges()[srcFacei];
        triEngine.triangulate
        (
            UIndirectList<point>(this->srcPatch_.points(), srcFacePoints)
        );
        srcTris.append(triEngine.triPoints(srcFacePoints));
        srcTriSrcFaces.resize(srcTris.size(), srcFacei);
        srcTriSrcEdges.append(triEngine.triEdges(srcFaceEdges));
    }
    forAll(this->tgtPatch_, tgtFacei)
    {
        const face tgtFacePoints = this->tgtPatch_[tgtFacei].reverseFace();
        const labelList tgtFaceEdges =
            reverseList(this->tgtPatch_.faceEdges()[tgtFacei]);
        triEngine.triangulate
        (
            UIndirectList<point>(this->tgtPatch_.points(), tgtFacePoints)
        );
        tgtTris.append(triEngine.triPoints(tgtFacePoints));
        tgtTriTgtFaces.resize(tgtTris.size(), tgtFacei);
        tgtTriTgtEdges.append(triEngine.triEdges(tgtFaceEdges));
        const label n = triEngine.triPoints().size();
        for (label i = 0; i < n; ++ i)
        {
            tgtTris[tgtTris.size() - n + i] =
                tgtTris[tgtTris.size() - n + i].reverseFace();
            tgtTriTgtEdges[tgtTris.size() - n + i] =
                reverseList(tgtTriTgtEdges[tgtTris.size() - n + i]);
        }
    }
    const primitiveTriPatch srcTriPatch
    (
        SubList<triFace>(srcTris, srcTris.size()),
        this->srcPatch_.points()
    );
    const primitiveTriPatch tgtTriPatch
    (
        SubList<triFace>(tgtTris, tgtTris.size()),
        this->tgtPatch_.points()
    );

    labelList srcTriPatchPointSrcPatchPoints(srcTriPatch.nPoints());
    labelList tgtTriPatchPointTgtPatchPoints(tgtTriPatch.nPoints());
    {
        labelList map
        (
            max(max(srcPatch.meshPoints()), max(tgtPatch.meshPoints())) + 1,
            -1
        );

        UIndirectList<label>(map, srcPatch.meshPoints()) =
            identityMap(srcPatch.nPoints());

        srcTriPatchPointSrcPatchPoints =
            renumber(map, srcTriPatch.meshPoints());

        UIndirectList<label>(map, srcPatch.meshPoints()) = -1;

        UIndirectList<label>(map, tgtPatch.meshPoints()) =
            identityMap(tgtPatch.nPoints());

        tgtTriPatchPointTgtPatchPoints =
            renumber(map, tgtTriPatch.meshPoints());

        UIndirectList<label>(map, tgtPatch.meshPoints()) = -1;
    }

    labelList srcTriPatchEdgeSrcPatchEdges(srcTriPatch.nEdges());
    labelList tgtTriPatchEdgeTgtPatchEdges(tgtTriPatch.nEdges());
    forAll(srcTriPatch, srcTrii)
    {
        forAll(srcTriPatch[srcTrii], srcTriEdgei)
        {
            const label srcTriPatchEdgei =
                srcTriPatch.faceEdges()[srcTrii][srcTriEdgei];
            const label srcPatchEdgei = srcTriSrcEdges[srcTrii][srcTriEdgei];
            srcTriPatchEdgeSrcPatchEdges[srcTriPatchEdgei] = srcPatchEdgei;
        }
    }
    forAll(tgtTriPatch, tgtTrii)
    {
        forAll(tgtTriPatch[tgtTrii], tgtTriEdgei)
        {
            const label tgtTriPatchEdgei =
                tgtTriPatch.faceEdges()[tgtTrii][tgtTriEdgei];
            const label tgtPatchEdgei = tgtTriTgtEdges[tgtTrii][tgtTriEdgei];
            tgtTriPatchEdgeTgtPatchEdges[tgtTriPatchEdgei] = tgtPatchEdgei;
        }
    }

    // Step 2: Do the tri-patch intersection
    const TriPatchIntersection<primitiveTriPatch, primitiveTriPatch> tpi
    (
        srcTriPatch,
        vectorField(srcPointNormals, srcTriPatchPointSrcPatchPoints),
        tgtTriPatch,
        snapTol
    );

    // Step 3: Map points and point associations back into the class data
    {
        this->points_ = tpi.points();

        // Point-point connectivity
        this->srcPointPoints_ =
            reorder
            (
                this->srcPatch().nPoints(),
                -1,
                srcTriPatchPointSrcPatchPoints,
                tpi.srcPointPoints()
            );
        this->tgtPointPoints_ =
            reorder
            (
                this->tgtPatch().nPoints(),
                -1,
                tgtTriPatchPointTgtPatchPoints,
                tpi.tgtPointPoints()
            );
        this->pointSrcPoints_ = tpi.pointSrcPoints();
        inplaceRenumber(srcTriPatchPointSrcPatchPoints, this->pointSrcPoints_);
        this->pointTgtPoints_ = tpi.pointTgtPoints();
        inplaceRenumber(tgtTriPatchPointTgtPatchPoints, this->pointTgtPoints_);

        // Point-edge connectivity
        this->srcEdgePoints_ =
            reorder
            (
                this->srcPatch_.nEdges(),
                DynamicList<label>(),
                srcTriPatchEdgeSrcPatchEdges,
                tpi.srcEdgePoints()
            );
        this->srcEdgePoints_.resize(this->srcPatch_.nEdges());
        this->tgtEdgePoints_ =
            reorder
            (
                this->tgtPatch_.nEdges(),
                DynamicList<label>(),
                tgtTriPatchEdgeTgtPatchEdges,
                tpi.tgtEdgePoints()
            );
        this->tgtEdgePoints_.resize(this->tgtPatch_.nEdges());
        this->pointSrcEdges_ = tpi.pointSrcEdges();
        inplaceRenumber(srcTriPatchEdgeSrcPatchEdges, this->pointSrcEdges_);
        this->pointTgtEdges_ = tpi.pointTgtEdges();
        inplaceRenumber(tgtTriPatchEdgeTgtPatchEdges, this->pointTgtEdges_);

        // Point-tri connectivity
        this->pointSrcFaces_ = tpi.pointSrcFaces();
        inplaceRenumber(srcTriSrcFaces, this->pointSrcFaces_);
        this->pointTgtFaces_ = tpi.pointTgtFaces();
        inplaceRenumber(tgtTriTgtFaces, this->pointTgtFaces_);

        // Point-face-diagonal connectivity
        auto setFaceDiagonalPoints = []
        (
            const List<DynamicList<label>>& edgePoints,
            const labelListList& edgeTris,
            const labelUList& triFaces,
            DynamicList<label>& pointFaces
        )
        {
            forAll(edgePoints, edgei)
            {
                const labelList& triis = edgeTris[edgei];

                if (triis.size() != 2)
                {
                    continue;
                }

                const label facei0 = triFaces[triis[0]];
                const label facei1 = triFaces[triis[1]];

                if (facei0 != facei1)
                {
                    continue;
                }

                const label facei = facei0;

                for
                (
                    label edgePointi = 1;
                    edgePointi < edgePoints[edgei].size() - 1;
                    ++ edgePointi
                )
                {
                    const label pointi = edgePoints[edgei][edgePointi];

                    pointFaces[pointi] = facei;
                }
            }
        };
        setFaceDiagonalPoints
        (
            tpi.srcEdgePoints(),
            srcTriPatch.edgeFaces(),
            srcTriSrcFaces,
            this->pointSrcFaces_
        );
        setFaceDiagonalPoints
        (
            tpi.tgtEdgePoints(),
            tgtTriPatch.edgeFaces(),
            tgtTriTgtFaces,
            this->pointTgtFaces_
        );
    }

    // Step 4: Construct faces and face associations by merging the faces that
    // result from the tri-patch intersection
    {
        // Loop the tri intersection faces, initialising a star at each one,
        // and propagating out to include everything connected with the same
        // source and target face. This is then part of the face intersection.
        boolList tpiFaceCombined(tpi.faces().size(), false);
        star tpiStar;
        forAll(tpi.faces(), tpiFacei)
        {
            // If this face has been done then skip to the next
            if (tpiFaceCombined[tpiFacei]) continue;

            // Get the source and target faces
            const label srcFacei =
                tpi.faceSrcFaces()[tpiFacei] != -1
              ? srcTriSrcFaces[tpi.faceSrcFaces()[tpiFacei]]
              : -1;
            const label tgtFacei =
                tpi.faceTgtFaces()[tpiFacei] != -1
              ? tgtTriTgtFaces[tpi.faceTgtFaces()[tpiFacei]]
              : -1;

            // Get the relevant set of face-edge and edge-face addressing
            const UList<labelList>& tpiFaceEdges = tpi.faceEdges();
            const UList<labelPair>& tpiEdgeFaces =
                srcFacei != -1 && tgtFacei != -1
              ? tpi.intersectEdgeFaces()
              : tpi.nonIntersectEdgeFaces();

            // Fill the star with all connected faces with the same source and
            // target face associations
            star::context tpiStarContext = tpiStar.populate
            (
                tpiFacei,
                true,
                [&](const label tpiEdgei, const label tpiFacei)
                {
                    const label tpiSrcFacei = tpi.faceSrcFaces()[tpiFacei];
                    const label tpiTgtFacei = tpi.faceTgtFaces()[tpiFacei];
                    return
                        (
                            tpiSrcFacei != -1
                          ? srcTriSrcFaces[tpiSrcFacei] == srcFacei
                          : srcFacei == -1
                        )
                     && (
                            tpiTgtFacei != -1
                          ? tgtTriTgtFaces[tpiTgtFacei] == tgtFacei
                          : tgtFacei == -1
                        );
                },
                tpiFaceEdges,
                tpiEdgeFaces
            );

            // Mark everything that has been visited
            forAllStarFaces(tpiStar, starFacei, tpiFacei)
            {
                tpiFaceCombined[tpiFacei] = true;
            }

            // Create the new face
            const label facei = this->faces_.size();
            this->faces_.append(face(tpiStar.starEdgeEdges().size()));
            if (srcFacei != -1)
            {
                this->srcFaceFaces_[srcFacei].append(facei);
            }
            if (tgtFacei != -1)
            {
                this->tgtFaceFaces_[tgtFacei].append(facei);
            }
            this->faceSrcFaces_.append(srcFacei);
            this->faceTgtFaces_.append(tgtFacei);
            forAllStarEdges(tpiStar, i, starEdgei, tpiEdgei)
            {
                const label tpiEdgeFacei =
                    tpiEdgeFaces[tpiEdgei][0] == -1
                 || tpiStar.faceStarFaces()[tpiEdgeFaces[tpiEdgei][0]] == -1;

                const label tpiFacei = tpiEdgeFaces[tpiEdgei][tpiEdgeFacei];
                const label tpiFaceEdgei =
                    findIndex(tpiFaceEdges[tpiFacei], tpiEdgei);

                this->faces_.last()[i] =
                    tpi.faces()[tpiFacei].faceEdge(tpiFaceEdgei).start();
            }
        }
    }

    // Step 5: Filter out points that aren't used
    {
        // Determine how many faces each point is a part of
        labelList pointNFaces(this->points_.size(), 0);
        forAll(this->faces_, facei)
        {
            const face& f = this->faces_[facei];
            forAll(f, fPointi)
            {
                ++ pointNFaces[f[fPointi]];
            }
        }

        // Remove points that are referenced by two faces or fewer, and which
        // do not correspond to original geometry or an intersection
        labelList oldPointNewPoints(this->points_.size(), -1);
        label pointi = 0;
        forAll(this->points_, pointj)
        {
            if
            (
                pointNFaces[pointj] > 2
             || (
                    pointNFaces[pointj] > 0
                 && (
                        this->pointSrcPoints_[pointj] != -1
                     || this->pointTgtPoints_[pointj] != -1
                     || (
                            this->pointSrcEdges_[pointj] != -1
                         && this->pointTgtEdges_[pointj] != -1
                        )
                    )
                )
            )
            {
                oldPointNewPoints[pointj] = pointi;
                this->points_[pointi] = this->points_[pointj];
                this->pointSrcPoints_[pointi] = this->pointSrcPoints_[pointj];
                this->pointTgtPoints_[pointi] = this->pointTgtPoints_[pointj];
                this->pointSrcEdges_[pointi] = this->pointSrcEdges_[pointj];
                this->pointTgtEdges_[pointi] = this->pointTgtEdges_[pointj];
                this->pointSrcFaces_[pointi] = this->pointSrcFaces_[pointj];
                this->pointTgtFaces_[pointi] = this->pointTgtFaces_[pointj];
                ++ pointi;
            }
        }
        this->points_.resize(pointi);
        this->pointSrcPoints_.resize(pointi);
        this->pointTgtPoints_.resize(pointi);
        this->pointSrcEdges_.resize(pointi);
        this->pointTgtEdges_.resize(pointi);
        this->pointSrcFaces_.resize(pointi);
        this->pointTgtFaces_.resize(pointi);

        // Map
        inplaceRenumber(oldPointNewPoints, this->faces_);
        inplaceRenumber(oldPointNewPoints, this->srcPointPoints_);
        inplaceRenumber(oldPointNewPoints, this->tgtPointPoints_);
        inplaceRenumber(oldPointNewPoints, this->srcEdgePoints_);
        inplaceRenumber(oldPointNewPoints, this->tgtEdgePoints_);

        // Remove deleted points
        inplaceRemove(this->faces_, label(-1));
        inplaceRemove(this->srcEdgePoints_, label(-1));
        inplaceRemove(this->tgtEdgePoints_, label(-1));
    }

    // Step 6: Reporting
    this->report();

    Info<< indent << this->faces_.size() << " couplings generated in "
            << time.cpuTimeIncrement() << 's' << endl;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SrcPatchType, class TgtPatchType>
Foam::FacePatchIntersection<SrcPatchType, TgtPatchType>::
~FacePatchIntersection()
{}


// ************************************************************************* //
