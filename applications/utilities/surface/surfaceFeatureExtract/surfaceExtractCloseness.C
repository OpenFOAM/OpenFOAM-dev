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

#include "surfaceFeatureExtract.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::extractCloseness
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

    const vectorField& normals = searchSurf.faceNormals();

    const scalar span = searchSurf.bounds().mag();

    const pointField start(searchSurf.faceCentres() - span*normals);
    const pointField end(searchSurf.faceCentres() + span*normals);
    const pointField& faceCentres = searchSurf.faceCentres();

    List<List<pointIndexHit>> allHitinfo;

    // Find all intersections (in order)
    searchSurf.findLineAll(start, end, allHitinfo);

    scalarField internalCloseness(start.size(), great);
    scalarField externalCloseness(start.size(), great);

    forAll(allHitinfo, fi)
    {
        const List<pointIndexHit>& hitInfo = allHitinfo[fi];

        processHit
        (
            internalCloseness[fi],
            externalCloseness[fi],
            fi,
            surf,
            start[fi],
            faceCentres[fi],
            end[fi],
            normals[fi],
            normals,
            hitInfo
        );
    }

    triSurfaceScalarField internalClosenessField
    (
        IOobject
        (
            sFeatFileName + ".internalCloseness",
            runTime.constant(),
            "triSurface",
            runTime
        ),
        surf,
        dimLength,
        internalCloseness
    );

    internalClosenessField.write();

    triSurfaceScalarField externalClosenessField
    (
        IOobject
        (
            sFeatFileName + ".externalCloseness",
            runTime.constant(),
            "triSurface",
            runTime
        ),
        surf,
        dimLength,
        externalCloseness
    );

    externalClosenessField.write();

    if (writeVTK)
    {
        const faceList faces(surf.faces());

        vtkSurfaceWriter().write
        (
            runTime.constantPath()/"triSurface",// outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "internalCloseness",                // fieldName
            internalCloseness,
            false,                              // isNodeValues
            true                                // verbose
        );

        vtkSurfaceWriter().write
        (
            runTime.constantPath()/"triSurface",// outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "externalCloseness",                // fieldName
            externalCloseness,
            false,                              // isNodeValues
            true                                // verbose
        );
    }
}


// ************************************************************************* //
