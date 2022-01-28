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

#include "zeroDimensionalFvMesh.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fvMesh Foam::zeroDimensionalFvMesh(const objectRegistry& db)
{
    pointField points(8);
    points[0] = vector(-0.5, -0.5, -0.5);
    points[1] = vector( 0.5, -0.5, -0.5);
    points[2] = vector( 0.5,  0.5, -0.5);
    points[3] = vector(-0.5,  0.5, -0.5);
    points[4] = vector(-0.5, -0.5,  0.5);
    points[5] = vector( 0.5, -0.5,  0.5);
    points[6] = vector( 0.5,  0.5,  0.5);
    points[7] = vector(-0.5,  0.5,  0.5);

    faceList faces = cellModeller::lookup("hex")->modelFaces();

    labelList owner(6, label(0));
    labelList neighbour(0);

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            db.time().timeName(),
            db,
            IOobject::READ_IF_PRESENT
        ),
        move(points),
        move(faces),
        move(owner),
        move(neighbour)
    );

    List<polyPatch*> patches
    (
        1,
        new emptyPolyPatch
        (
            "boundary",
            6,
            0,
            0,
            mesh.boundaryMesh(),
            emptyPolyPatch::typeName
        )
    );

    mesh.addFvPatches(patches);

    return mesh;
}

// ************************************************************************* //
