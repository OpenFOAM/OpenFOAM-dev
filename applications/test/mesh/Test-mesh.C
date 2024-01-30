/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh, no clear-out\n" << endl;
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.name(),
            runTime,
            IOobject::MUST_READ
        ),
        false
    );

    Info<< mesh.C() << endl;
    Info<< mesh.V() << endl;

    surfaceVectorField Cf = mesh.Cf();

    Info<< Cf << endl;

    // Test construct from cellShapes
    {
        pointField points(mesh.points());
        cellShapeList shapes(mesh.cellShapes());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        faceListList boundaryFaces(pbm.size());
        forAll(pbm, patchi)
        {
            boundaryFaces[patchi] = pbm[patchi];
        }
        wordList boundaryPatchNames(pbm.names());
        PtrList<dictionary> boundaryDicts(pbm.size());
        forAll(pbm, patchi)
        {
            OStringStream os;
            os << pbm[patchi];
            IStringStream is(os.str());
            boundaryDicts.set(patchi, new dictionary(is));
        }

        word defaultBoundaryPatchName = "defaultFaces";
        word defaultBoundaryPatchType = emptyPolyPatch::typeName;

        fvMesh newMesh
        (
            IOobject
            (
                "newMesh",
                runTime.name(),
                runTime,
                Foam::IOobject::NO_READ
            ),
            clone(points),
            shapes,
            boundaryFaces,
            boundaryPatchNames,
            boundaryDicts,
            defaultBoundaryPatchName,
            defaultBoundaryPatchType
        );

        Info<< newMesh.C() << endl;
        Info<< newMesh.V() << endl;

        surfaceVectorField Cf = newMesh.Cf();

        Info<< Cf << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
