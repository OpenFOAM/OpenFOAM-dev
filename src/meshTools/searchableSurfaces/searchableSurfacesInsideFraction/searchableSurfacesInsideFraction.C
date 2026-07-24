/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "searchableSurfacesInsideFraction.H"
#include "syncTools.H"
#include "cutPolyIntegral.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::searchableSurfaces::insideFraction
(
    const searchableSurface& surface,
    const polyMesh& mesh
)
{
    const boundBox& surfaceBb = surface.bounds();

    // Step 1. Calculate a length scale for each mesh point. This is used to
    // determine if a point is in the range of the surface and thereby prevent
    // expensive surface tests on points that are far away.

    scalarField pointLengthScale(mesh.nPoints(), -vGreat);
    forAll(mesh.faces(), facei)
    {
        forAll(mesh.faces()[facei], faceEdgei)
        {
            const edge e = mesh.faces()[facei].faceEdge(faceEdgei);
            const scalar l = e.mag(mesh.points());

            forAll(e, i)
            {
                pointLengthScale[e[i]] = max(pointLengthScale[e[i]], l);
                pointLengthScale[e[i]] = max(pointLengthScale[e[i]], l);
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        pointLengthScale,
        maxEqOp(),
        -vGreat
    );

    // Step 2. Calculate a signed distance field from the surface

    // Determine points in range of the surface
    DynamicList<label> nearPointPointis(mesh.nPoints());
    forAll(mesh.points(), pointi)
    {
        const point& p = mesh.points()[pointi];
        const vector l = pointLengthScale[pointi]*vector::one;

        if (surfaceBb.overlaps(boundBox(p - l, p + l)))
        {
            nearPointPointis.append(pointi);
        }
    }
    nearPointPointis.shrink();
    const pointField nearPoints(mesh.points(), nearPointPointis);

    // Search for mesh points' nearest surface points
    List<pointIndexHit> nearPointHits(nearPoints.size());
    surface.findNearest
    (
        nearPoints,
        scalarField(nearPoints.size(), rootVGreat),
        nearPointHits
    );

    // Determine whether the mesh points are inside or outside the surface
    List<volumeType> nearPointsVolumeType(nearPoints.size());
    surface.getVolumeType(nearPoints, nearPointsVolumeType);

    // Convert points into distances, and use the inside/outside status to
    // assign a sign to that distance
    scalarField pointDistances(mesh.nPoints(), rootVGreat);
    forAll(nearPoints, nearPointi)
    {
        const label pointi = nearPointPointis[nearPointi];

        const scalar d =
            mag(nearPoints[nearPointi] - nearPointHits[nearPointi].hitPoint());

        switch (nearPointsVolumeType[nearPointi])
        {
            case volumeType::unknown:
            case volumeType::mixed:
                pointDistances[pointi] = vGreat;
                break;
            case volumeType::outside:
                pointDistances[pointi] = d;
                break;
            case volumeType::inside:
                pointDistances[pointi] = -d;
                break;
        }
    }

    // Step 3. Cut the mesh at the zero iso-surface of the signed distance
    // field, and obtain the volumes that result.

    const cellEdgeAddressingList& cAddrs = cellEdgeAddressingList::New(mesh);

    List<List<labelPair>> faceCuts(mesh.faces().size());
    forAll(mesh.faces(), facei)
    {
        faceCuts[facei] =
            cutPoly::faceCuts
            (
                mesh.faces()[facei],
                pointDistances,
                0
            );
    }

    List<labelListList> cellCuts(mesh.cells().size());
    forAll(mesh.cells(), celli)
    {
        cellCuts[celli] =
            cutPoly::cellCuts
            (
                mesh.cells()[celli],
                cAddrs[celli],
                mesh.faces(),
                faceCuts,
                pointDistances,
                0
            );
    }

    vectorField faceCutAreas(mesh.faces().size());
    pointField faceCutCentres(mesh.faces().size());
    forAll(mesh.faces(), facei)
    {
        const Tuple2<vector, tensor> fCutAreaCentre =
            cutPoly::faceCutAreaIntegral
            (
                mesh.faces()[facei],
                mesh.faceAreas()[facei],
                mesh.faceCentres()[facei],
                faceCuts[facei],
                mesh.points(),
                mesh.points(),
                pointDistances,
                0,
                true
            );

        const vector fCutArea = fCutAreaCentre.first();
        const scalar fMagSqrCutArea = magSqr(fCutArea);

        faceCutAreas[facei] = fCutArea;
        faceCutCentres[facei] =
            fMagSqrCutArea > vSmall
          ? (fCutArea & fCutAreaCentre.second())/fMagSqrCutArea
          : mesh.faceCentres()[facei];
    }

    scalarField cellCutVolumes(mesh.cells().size());
    forAll(mesh.cells(), celli)
    {
        cellCutVolumes[celli] =
            cutPoly::cellCutVolume
            (
                mesh.cells()[celli],
                cAddrs[celli],
                mesh.cellVolumes()[celli],
                cellCuts[celli],
                mesh.faces(),
                mesh.faceAreas(),
                mesh.faceCentres(),
                faceCutAreas,
                mesh.points(),
                pointDistances,
                0,
                true
            );
    }

    // Step 4. Divide the cut volumes through by the total cell volumes to get
    // the volume fraction

    return cellCutVolumes/mesh.cellVolumes();
}


Foam::tmp<Foam::scalarField> Foam::searchableSurfaces::insideFraction
(
    const searchableSurface& surface,
    const polyPatch& patch
)
{
    const pointField& points = patch.localPoints();
    const faceList& faces = patch.localFaces();
    const vectorField::subField& faceAreas = patch.faceAreas();

    // Search for mesh points' nearest surface points
    List<pointIndexHit> pointHits(points.size());
    surface.findNearest
    (
        points,
        scalarField(points.size(), rootVGreat),
        pointHits
    );

    // Determine whether the mesh points are inside or outside the surface
    List<volumeType> pointsVolumeType(points.size());
    surface.getVolumeType(points, pointsVolumeType);

    // Convert points into distances, and use the inside/outside status to
    // assign a sign to that distance
    scalarField pointDistances(points.size(), rootVGreat);
    forAll(points, pointi)
    {
        const scalar d =
            mag(points[pointi] - pointHits[pointi].hitPoint());

        switch (pointsVolumeType[pointi])
        {
            case volumeType::unknown:
            case volumeType::mixed:
                pointDistances[pointi] = vGreat;
                break;
            case volumeType::outside:
                pointDistances[pointi] = d;
                break;
            case volumeType::inside:
                pointDistances[pointi] = -d;
                break;
        }
    }

    // Cut the faces at the zero iso-surface of the signed distance field, and
    // obtain the areas that result.
    List<List<labelPair>> faceCuts(faces.size());
    forAll(faces, facei)
    {
        faceCuts[facei] =
            cutPoly::faceCuts
            (
                faces[facei],
                pointDistances,
                0
            );
    }

    vectorField faceCutAreas(faces.size());
    forAll(faces, facei)
    {
        faceCutAreas[facei] =
            cutPoly::faceCutArea
            (
                faces[facei],
                faceAreas[facei],
                faceCuts[facei],
                points,
                pointDistances,
                0,
                true
            );
    }

    // Divide the cut areas through by the total face areas to get the area
    // fraction
    return mag(faceCutAreas)/mag(faceAreas);
}


// ************************************************************************* //
