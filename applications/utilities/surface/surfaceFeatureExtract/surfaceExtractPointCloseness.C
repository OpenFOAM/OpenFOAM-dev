/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "surfaceFeatureExtract.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::processHit
(
    scalar& internalCloseness,
    scalar& externalCloseness,
    const label fi,
    const triSurface& surf,
    const point& start,
    const point& p,
    const point& end,
    const vector& normal,
    const vectorField& normals,
    const List<pointIndexHit>& hitInfo
)
{
    if (hitInfo.size() < 1)
    {
        drawHitProblem(fi, surf, start, p, end, hitInfo);
    }
    else if (hitInfo.size() == 1)
    {
        if (!hitInfo[0].hit())
        {
        }
        else if (hitInfo[0].index() != fi)
        {
            drawHitProblem(fi, surf, start, p, end, hitInfo);
        }
    }
    else
    {
        label ownHiti = -1;

        forAll(hitInfo, hI)
        {
            // Find the hit on the triangle that launched the ray

            if (hitInfo[hI].index() == fi)
            {
                ownHiti = hI;
                break;
            }
        }

        if (ownHiti < 0)
        {
            drawHitProblem(fi, surf, start, p, end, hitInfo);
        }
        else if (ownHiti == 0)
        {
            // There are no internal hits, the first hit is the
            // closest external hit

            if
            (
                (normal & normals[hitInfo[ownHiti + 1].index()])
              < externalToleranceCosAngle
            )
            {
                externalCloseness = min
                (
                    externalCloseness,
                    mag(p - hitInfo[ownHiti + 1].hitPoint())
                );
            }
        }
        else if (ownHiti == hitInfo.size() - 1)
        {
            // There are no external hits, the last but one hit is
            // the closest internal hit

            if
            (
                (normal & normals[hitInfo[ownHiti - 1].index()])
              < internalToleranceCosAngle
            )
            {
                internalCloseness = min
                (
                    internalCloseness,
                    mag(p - hitInfo[ownHiti - 1].hitPoint())
                );
            }
        }
        else
        {
            if
            (
                (normal & normals[hitInfo[ownHiti + 1].index()])
              < externalToleranceCosAngle
            )
            {
                externalCloseness = min
                (
                    externalCloseness,
                    mag(p - hitInfo[ownHiti + 1].hitPoint())
                );
            }

            if
            (
                (normal & normals[hitInfo[ownHiti - 1].index()])
              < internalToleranceCosAngle
            )
            {
                internalCloseness = min
                (
                    internalCloseness,
                    mag(p - hitInfo[ownHiti - 1].hitPoint())
                );
            }
        }
    }
}


void Foam::extractPointCloseness
(
    const fileName &sFeatFileName,
    const Time& runTime,
    const triSurface &surf,
    const bool writeVTK
)
{
    // Searchable triSurface
    const triSurfaceMesh searchSurf
    (
        IOobject
        (
            sFeatFileName + ".closeness",
            runTime.constant(),
            "triSurface",
            runTime
        ),
        surf
    );


    // Prepare start and end points for intersection tests

    const pointField& points = searchSurf.points();
    const labelList& meshPoints = searchSurf.meshPoints();
    const pointField& faceCentres = searchSurf.faceCentres();
    const vectorField& normals = searchSurf.faceNormals();
    const labelListList& pointFaces = searchSurf.pointFaces();

    const scalar span = searchSurf.bounds().mag();

    label nPointFaces = 0;
    forAll(pointFaces, pfi)
    {
        nPointFaces += pointFaces[pfi].size();
    }

    pointField facePoints(nPointFaces);
    pointField start(nPointFaces);
    pointField end(nPointFaces);

    label i = 0;
    forAll(points, pi)
    {
        forAll(pointFaces[pi], pfi)
        {
            const label fi = pointFaces[pi][pfi];

            facePoints[i] = (0.9*points[meshPoints[pi]] + 0.1*faceCentres[fi]);
            const vector& n = normals[fi];

            start[i] = facePoints[i] - span*n;
            end[i] = facePoints[i] + span*n;

            i++;
        }
    }

    List<List<pointIndexHit>> allHitinfo;

    // Find all intersections (in order)
    searchSurf.findLineAll(start, end, allHitinfo);

    scalarField internalCloseness(points.size(), great);
    scalarField externalCloseness(points.size(), great);

    i = 0;
    forAll(points, pi)
    {
        forAll(pointFaces[pi], pfi)
        {
            const label fi = pointFaces[pi][pfi];
            const List<pointIndexHit>& hitInfo = allHitinfo[i];

            processHit
            (
                internalCloseness[pi],
                externalCloseness[pi],
                fi,
                surf,
                start[i],
                facePoints[i],
                end[i],
                normals[fi],
                normals,
                hitInfo
            );

            i++;
        }
    }

    triSurfacePointScalarField internalClosenessPointField
    (
        IOobject
        (
            sFeatFileName + ".internalPointCloseness",
            runTime.constant(),
            "triSurface",
            runTime
        ),
        surf,
        dimLength,
        internalCloseness
    );

    internalClosenessPointField.write();

    triSurfacePointScalarField externalClosenessPointField
    (
        IOobject
        (
            sFeatFileName + ".externalPointCloseness",
            runTime.constant(),
            "triSurface",
            runTime
        ),
        surf,
        dimLength,
        externalCloseness
    );

    externalClosenessPointField.write();

    if (writeVTK)
    {
        const faceList faces(surf.faces());
        const Map<label>& meshPointMap = surf.meshPointMap();

        forAll(meshPointMap, pi)
        {
            internalCloseness[pi] =
                internalClosenessPointField[meshPointMap[pi]];

            externalCloseness[pi] =
                externalClosenessPointField[meshPointMap[pi]];
        }

        vtkSurfaceWriter().write
        (
            runTime.constantPath()/"triSurface",// outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "internalPointCloseness",           // fieldName
            internalCloseness,
            true,                               // isNodeValues
            true                                // verbose
        );

        vtkSurfaceWriter().write
        (
            runTime.constantPath()/"triSurface",// outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "externalPointCloseness",           // fieldName
            externalCloseness,
            true,                               // isNodeValues
            true                                // verbose
        );
    }
}


// ************************************************************************* //
