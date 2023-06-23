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

Application
    Test-mappedPatch

Description
    Test mapped patches and boundary conditions by mapping cell and face
    centres and writing out the resulting point connections

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "meshTools.H"
#include "Time.H"
#include "OBJstream.H"
#include "volFields.H"
#include "mappedPolyPatch.H"
#include "mappedInternalPolyPatch.H"
#include "mappedValueFvPatchFields.H"
#include "mappedInternalValueFvPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    timeSelector::select0(runTime, args);
    #include "createMesh.H"

    wordList patchFieldTypes
    (
        mesh.boundaryMesh().size(),
        calculatedFvPatchVectorField::typeName
    );

    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isA<mappedPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            patchFieldTypes[patchi] =
                mappedValueFvPatchVectorField::typeName;
        }
        if (isA<mappedInternalPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            patchFieldTypes[patchi] =
                mappedInternalValueFvPatchVectorField::typeName;
        }
    }

    volVectorField cc
    (
        IOobject
        (
            "cc",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C(),
        patchFieldTypes
    );

    cc.correctBoundaryConditions();

    forAll(cc.boundaryField(), patchi)
    {
        const fvPatchVectorField& ccp = cc.boundaryField()[patchi];

        if
        (
            isA<mappedValueFvPatchVectorField>(ccp)
         || isA<mappedInternalValueFvPatchVectorField>(ccp)
        )
        {
            OBJstream obj(mesh.boundaryMesh()[patchi].name() + ".obj");

            Pout<< "Detected a " << ccp.patch().type() << " field on patch \""
                << ccp.patch().name() << "\". Writing point connections to "
                << obj.name() << "." << endl;

            forAll(ccp, i)
            {
                obj.write(linePointRef(ccp.patch().Cf()[i], ccp[i]));
            }
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
