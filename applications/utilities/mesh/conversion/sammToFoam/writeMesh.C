/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Create intermediate mesh files from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sammMesh::writeMesh()
{
    if (isShapeMesh_)
    {
        Info<< "This is a shapeMesh." << endl;

        polyMesh pShapeMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_
            ),
            xferCopy(points_),           // we could probably re-use the data
            cellShapes_,
            boundary_,
            patchNames_,
            patchTypes_,
            defaultFacesName_,
            defaultFacesType_,
            patchPhysicalTypes_
        );

        Info<< "Writing polyMesh" << endl;
        pShapeMesh.write();
    }
    else
    {
        // This is a polyMesh.

        createPolyMeshData();

        Info<< "This is a polyMesh" << endl;

        polyMesh pMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_
            ),
            xferCopy(points_),           // we could probably re-use the data
            xferCopy(meshFaces_),
            xferCopy(cellPolys_)
        );

        pMesh.addPatches(polyBoundaryPatches(pMesh));

        Info<< "Writing polyMesh" << endl;
        pMesh.write();
    }
}


// ************************************************************************* //
