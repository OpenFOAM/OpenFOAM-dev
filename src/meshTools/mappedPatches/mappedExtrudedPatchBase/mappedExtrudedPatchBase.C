/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "mappedExtrudedPatchBase.H"
#include "LayerInfoData.H"
#include "PointEdgeLayerInfoData.H"
#include "FaceCellWave.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedExtrudedPatchBase, 0);

    template<>
    inline bool contiguous<LayerInfoData<Pair<vector>>>()
    {
        return true;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::mappedExtrudedPatchBase::patchFaceAreas() const
{
    if (isExtrudedRegion_)
    {
        if (!bottomFaceAreasPtr_.valid())
        {
            const polyMesh& mesh = patch_.boundaryMesh().mesh();
            const polyPatch& pp = patch_;

            // If this is the extruded region we need to work out what the
            // corresponding areas and centres are on the bottom patch. We do
            // this by waving these values across the layers.

            // Initialise layer data on the patch faces
            labelList initialFaces1(pp.size());
            List<layerInfo> initialFaceInfo1(pp.size());
            forAll(pp, ppFacei)
            {
                initialFaces1[ppFacei] = pp.start() + ppFacei;
                initialFaceInfo1[ppFacei] = layerInfo(0, -1);
            }

            // Wave across the mesh layers
            List<layerInfo> faceInfo1(mesh.nFaces());
            List<layerInfo> cellInfo1(mesh.nCells());
            FaceCellWave<layerInfo> wave1(mesh, faceInfo1, cellInfo1);
            wave1.setFaceInfo(initialFaces1, initialFaceInfo1);
            const label nIterations1 =
                wave1.iterate(mesh.globalData().nTotalCells() + 1);

            // Count how many opposite faces the wave ended on
            label nInitialFaces2 = 0;
            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                ++ facei
            )
            {
                if
                (
                    faceInfo1[facei].valid(wave1.data())
                 && faceInfo1[facei].faceLayer() == nIterations1
                )
                {
                    nInitialFaces2 ++;
                }
            }

            // Initialise data on the opposite faces. Store the area and centre.
            labelList initialFaces2(nInitialFaces2);
            List<LayerInfoData<Pair<vector>>> initialFaceInfo2(nInitialFaces2);
            label initialFace2i = 0;
            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                ++ facei
            )
            {
                if
                (
                    faceInfo1[facei].valid(wave1.data())
                 && faceInfo1[facei].faceLayer() != 0
                )
                {
                    initialFaces2[initialFace2i] = facei;
                    initialFaceInfo2[initialFace2i] =
                        LayerInfoData<Pair<vector>>
                        (
                            0,
                            -1,
                            Pair<vector>
                            (
                                mesh.faceAreas()[facei],
                                mesh.faceCentres()[facei]
                            )
                        );
                    initialFace2i ++;
                }
            }

            // Wave back across the mesh layers
            List<LayerInfoData<Pair<vector>>> faceInfo2(mesh.nFaces());
            List<LayerInfoData<Pair<vector>>> cellInfo2(mesh.nCells());
            FaceCellWave<LayerInfoData<Pair<vector>>> wave2
            (
                mesh,
                initialFaces2,
                initialFaceInfo2,
                faceInfo2,
                cellInfo2,
                mesh.globalData().nTotalCells() + 1
            );

            // Unpack into this patch's bottom face areas and centres. Note
            // that the face area needs flipping as it relates to a patch on
            // the other side of the extruded region.
            bottomFaceAreasPtr_.set(new vectorField(pp.size()));
            bottomFaceCentresPtr_.set(new pointField(pp.size()));
            forAll(pp, ppFacei)
            {
                const LayerInfoData<Pair<vector>>& info =
                    faceInfo2[pp.start() + ppFacei];

                static nil td;

                if (!info.valid(td))
                {
                    FatalErrorInFunction
                        << "Mesh \"" << mesh.name()
                        << "\" is not layered from the extruded patch "
                        << "\"" << pp.name() << "\"" << exit(FatalError);
                }

                bottomFaceAreasPtr_()[ppFacei] = - info.data().first();
                bottomFaceCentresPtr_()[ppFacei] = info.data().second();
            }
        }

        return bottomFaceAreasPtr_();
    }

    return mappedPatchBase::patchFaceAreas();
}


Foam::tmp<Foam::pointField>
Foam::mappedExtrudedPatchBase::patchFaceCentres() const
{
    if (isExtrudedRegion_)
    {
        if (!bottomFaceCentresPtr_.valid())
        {
            patchFaceAreas();
        }

        return bottomFaceCentresPtr_();
    }

    return mappedPatchBase::patchFaceCentres();
}


Foam::tmp<Foam::pointField>
Foam::mappedExtrudedPatchBase::patchLocalPoints() const
{
    if (isExtrudedRegion_)
    {
        if (!bottomLocalPointsPtr_.valid())
        {
            const polyMesh& mesh = patch_.boundaryMesh().mesh();
            const polyPatch& pp = patch_;

            // If this is the extruded region we need to work out what the
            // corresponding points are on the bottom patch. We do this by
            // waving these location across the layers.

            // Initialise layer data on the patch points
            labelList initialPoints1(pp.nPoints());
            List<pointEdgeLayerInfo> initialPointInfo1(pp.nPoints());
            forAll(pp.meshPoints(), ppPointi)
            {
                initialPoints1[ppPointi] = pp.meshPoints()[ppPointi];
                initialPointInfo1[ppPointi] = pointEdgeLayerInfo(0);
            }

            // Wave across the mesh layers
            List<pointEdgeLayerInfo> pointInfo1(mesh.nPoints());
            List<pointEdgeLayerInfo> edgeInfo1(mesh.nEdges());
            PointEdgeWave<pointEdgeLayerInfo> wave1
            (
                mesh,
                pointInfo1,
                edgeInfo1
            );
            wave1.setPointInfo(initialPoints1, initialPointInfo1);
            const label nIterations1 =
                wave1.iterate(mesh.globalData().nTotalPoints() + 1);

            if (debug)
            {
                pointScalarField pointLayer
                (
                    pointScalarField::New
                    (
                        typedName("pointLayer"),
                        pointMesh::New(mesh),
                        dimensionedScalar(dimless, 0)
                    )
                );
                forAll(pointInfo1, pointi)
                {
                    pointLayer[pointi] =
                        pointInfo1[pointi].valid(wave1.data())
                      ? pointInfo1[pointi].pointLayer()
                      : -1;
                }
                pointLayer.write();
            }

            // Count how many opposite points the wave ended on
            label nInitialPoints2 = 0;
            forAll(pointInfo1, pointi)
            {
                if
                (
                    pointInfo1[pointi].valid(wave1.data())
                 && pointInfo1[pointi].pointLayer() == nIterations1
                )
                {
                    nInitialPoints2 ++;
                }
            }

            // Initialise data on the opposite points. Store the position.
            labelList initialPoints2(nInitialPoints2);
            List<PointEdgeLayerInfoData<point>>
                initialPointInfo2(nInitialPoints2);
            label initialPoint2i = 0;
            forAll(pointInfo1, pointi)
            {
                if
                (
                    pointInfo1[pointi].valid(wave1.data())
                 && pointInfo1[pointi].pointLayer() == nIterations1
                )
                {
                    initialPoints2[initialPoint2i] = pointi;
                    initialPointInfo2[initialPoint2i] =
                        PointEdgeLayerInfoData<point>(0, mesh.points()[pointi]);
                    initialPoint2i ++;
                }
            }

            // Wave back across the mesh layers
            List<PointEdgeLayerInfoData<point>> pointInfo2(mesh.nPoints());
            List<PointEdgeLayerInfoData<point>> edgeInfo2(mesh.nEdges());
            PointEdgeWave<PointEdgeLayerInfoData<point>> wave2
            (
                mesh,
                initialPoints2,
                initialPointInfo2,
                pointInfo2,
                edgeInfo2,
                mesh.globalData().nTotalCells() + 1
            );

            // Unpack into this patch's bottom local points
            bottomLocalPointsPtr_.set(new pointField(pp.nPoints()));
            forAll(pp.meshPoints(), ppPointi)
            {
                const PointEdgeLayerInfoData<point>& info =
                    pointInfo2[pp.meshPoints()[ppPointi]];

                static nil td;

                if (!info.valid(td))
                {
                    FatalErrorInFunction
                        << "Mesh \"" << mesh.name()
                        << "\" is not layered from the extruded patch "
                        << "\"" << pp.name() << "\"" << exit(FatalError);
                }

                bottomLocalPointsPtr_()[ppPointi] = info.data();
            }

            if (debug)
            {
                pointVectorField pointOffset
                (
                    pointVectorField::New
                    (
                        typedName("pointOffset"),
                        pointMesh::New(mesh),
                        dimensionedVector(dimLength, Zero)
                    )
                );
                forAll(pp.meshPoints(), ppPointi)
                {
                    pointOffset[pp.meshPoints()[ppPointi]] =
                        bottomLocalPointsPtr_()[ppPointi]
                      - pp.localPoints()[ppPointi];
                }
                pointOffset.write();
            }
        }

        return bottomLocalPointsPtr_();
    }

    return mappedPatchBase::patchLocalPoints();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase(const polyPatch& pp)
:
    mappedPatchBase(pp),
    isExtrudedRegion_(false)
{}


Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase
(
    const polyPatch& pp,
    const word& nbrRegionName,
    const word& nbrPatchName,
    const bool isExtrudedRegion,
    const cyclicTransform& transform
)
:
    mappedPatchBase(pp, nbrRegionName, nbrPatchName, transform),
    isExtrudedRegion_(isExtrudedRegion)
{}


Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict,
    const bool defaultTransformIsNone
)
:
    mappedPatchBase(pp, dict, defaultTransformIsNone),
    isExtrudedRegion_(dict.lookup<bool>("isExtrudedRegion"))
{}


Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase
(
    const polyPatch& pp,
    const mappedExtrudedPatchBase& mepb
)
:
    mappedPatchBase(pp, mepb),
    isExtrudedRegion_(mepb.isExtrudedRegion_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedExtrudedPatchBase::~mappedExtrudedPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedExtrudedPatchBase::clearOut()
{
    mappedPatchBase::clearOut();
    if (reMapAfterMove_)
    {
        bottomFaceCentresPtr_.clear();
    }
}


void Foam::mappedExtrudedPatchBase::write(Ostream& os) const
{
    mappedPatchBase::write(os);
    writeEntry(os, "isExtrudedRegion", isExtrudedRegion_);
}


// ************************************************************************* //
