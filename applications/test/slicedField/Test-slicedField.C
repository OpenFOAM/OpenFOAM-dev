/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

Application
    Test-slicedField

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "SlicedGeometricField.H"
#include "slicedFvPatchFields.H"
#include "slicedSurfaceFields.H"
#include "surfaceInterpolate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Info<< min(p, p);

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    SlicedGeometricField<vector, volMesh> C
    (
        IOobject
        (
            "C2",
            runTime.name(),
            mesh
        ),
        mesh,
        dimLength,
        mesh.cellCentres(),
        mesh.faceCentres()
    );

    Info<< C << endl;
    Info<< (C & U) << endl;

    SlicedGeometricField<vector, surfaceMesh> Sf
    (
        IOobject
        (
            "Sf2",
            runTime.name(),
            mesh
        ),
        mesh,
        dimArea,
        mesh.faceAreas()
    );

    // Info<< Sf << endl;

    surfaceScalarField phi(fvc::interpolate(U) & mesh.Sf());

    scalarField completePhi
    (
        slicedSurfaceScalarField
        (
            IOobject
            (
                "slicedPhi",
                runTime.name(),
                mesh
            ),
            phi,
            false
        ).splice()
    );

    // Info<< completePhi << endl;

    return 0;
}


// ************************************************************************* //
