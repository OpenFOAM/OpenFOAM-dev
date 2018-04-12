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
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::triSurfacePointScalarField>>
Foam::extractPointCloseness
(
    const triSurfaceMesh& surf
)
{
    const Time& runTime = surf.objectRegistry::time();

    // Prepare start and end points for intersection tests

    const pointField& points = surf.points();
    const labelList& meshPoints = surf.meshPoints();
    const pointField& faceCentres = surf.faceCentres();
    const vectorField& normals = surf.faceNormals();
    const labelListList& pointFaces = surf.pointFaces();

    const scalar span = surf.bounds().mag();

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

    List<pointIndexHitList> allHitinfo;

    // Find all intersections (in order)
    surf.findLineAll(start, end, allHitinfo);

    scalarField internalCloseness(points.size(), great);
    scalarField externalCloseness(points.size(), great);

    i = 0;
    forAll(points, pi)
    {
        forAll(pointFaces[pi], pfi)
        {
            const label fi = pointFaces[pi][pfi];
            const pointIndexHitList& hitInfo = allHitinfo[i];

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

    return Pair<tmp<triSurfacePointScalarField>>
    (
        tmp<triSurfacePointScalarField>
        (
            new triSurfacePointScalarField
            (
                IOobject
                (
                    surf.objectRegistry::name() + ".internalPointCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf,
                dimLength,
                internalCloseness
            )
        ),

        tmp<triSurfacePointScalarField>
        (
            new triSurfacePointScalarField
            (
                IOobject
                (
                    surf.objectRegistry::name() + ".externalPointCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf,
                dimLength,
                externalCloseness
            )
        )
    );
}


// ************************************************************************* //
