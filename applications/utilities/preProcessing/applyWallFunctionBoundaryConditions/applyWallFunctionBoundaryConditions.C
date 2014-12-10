/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    applyWallFunctionBounaryConditions

Description
    Updates OpenFOAM RAS cases to use the new (v1.6) wall function framework.

    Attempts to determine whether case is compressible or incompressible, or
    can be supplied with -compressible command line argument.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "wallPolyPatch.H"

#include "incompressible/RAS/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.H"
#include "incompressible/RAS/derivedFvPatchFields/wallFunctions/kqRWallFunctions/kqRWallFunction/kqRWallFunctionFvPatchField.H"
#include "incompressible/RAS/derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkWallFunction/nutkWallFunctionFvPatchScalarField.H"
#include "incompressible/RAS/derivedFvPatchFields/wallFunctions/omegaWallFunctions/omegaWallFunction/omegaWallFunctionFvPatchScalarField.H"

#include "compressible/RAS/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.H"
#include "compressible/RAS/derivedFvPatchFields/wallFunctions/kqRWallFunctions/kqRWallFunction/kqRWallFunctionFvPatchField.H"
#include "compressible/RAS/derivedFvPatchFields/wallFunctions/mutWallFunctions/mutkWallFunction/mutkWallFunctionFvPatchScalarField.H"
#include "compressible/RAS/derivedFvPatchFields/wallFunctions/omegaWallFunctions/omegaWallFunction/omegaWallFunctionFvPatchScalarField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool caseIsCompressible(const fvMesh& mesh)
{
    // Attempt flux field
    IOobject phiHeader
    (
        "phi",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (phiHeader.headerOk())
    {
        surfaceScalarField phi(phiHeader, mesh);
        if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            return true;
        }
    }

    // Attempt density field
    IOobject rhoHeader
    (
        "rho",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (rhoHeader.headerOk())
    {
        volScalarField rho(rhoHeader, mesh);
        if (rho.dimensions() == dimDensity)
        {
            return true;
        }
    }

    // Attempt pressure field
    IOobject pHeader
    (
        "p",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (pHeader.headerOk())
    {
        volScalarField p(pHeader, mesh);
        if (p.dimensions() == dimMass/sqr(dimTime)/dimLength)
        {
            return true;
        }
    }

    // If none of the above are true, assume that the case is incompressible
    return false;
}


void createVolScalarField
(
    const fvMesh& mesh,
    const word& fieldName,
    const dimensionSet& dims
)
{
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!fieldHeader.headerOk())
    {
        Info<< "Creating field " << fieldName << nl << endl;

        volScalarField field
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dims, 0.0)
        );

        field.write();
    }
}


void replaceBoundaryType
(
    const fvMesh& mesh,
    const word& fieldName,
    const word& boundaryType,
    const string& boundaryValue
)
{
    IOobject header
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!header.headerOk())
    {
        return;
    }

    Info<< "Updating boundary types for field " << header.name() << endl;

    const word oldTypeName = IOdictionary::typeName;
    const_cast<word&>(IOdictionary::typeName) = word::null;

    IOdictionary dict(header);

    const_cast<word&>(IOdictionary::typeName) = oldTypeName;
    const_cast<word&>(dict.type()) = dict.headerClassName();

    // Make a backup of the old file
    if (mvBak(dict.objectPath(), "old"))
    {
        Info<< "    Backup original file to "
            << (dict.objectPath() + ".old") << endl;
    }

    // Loop through boundary patches and update
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    dictionary& boundaryDict = dict.subDict("boundaryField");
    forAll(bMesh, patchI)
    {
        if (isA<wallPolyPatch>(bMesh[patchI]))
        {
            word patchName = bMesh[patchI].name();
            dictionary& oldPatch = boundaryDict.subDict(patchName);

            dictionary newPatch(dictionary::null);
            newPatch.add("type", boundaryType);
            newPatch.add("value", ("uniform " + boundaryValue).c_str());

            oldPatch = newPatch;
        }
    }

    Info<< "    writing updated " << dict.name() << nl << endl;
    dict.regIOobject::write();
}


void updateCompressibleCase(const fvMesh& mesh)
{
    Info<< "Case treated as compressible" << nl << endl;
    createVolScalarField
    (
        mesh,
        "mut",
        dimArea/dimTime*dimDensity
    );
    replaceBoundaryType
    (
        mesh,
        "mut",
        compressible::mutkWallFunctionFvPatchScalarField::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "epsilon",
        compressible::epsilonWallFunctionFvPatchScalarField::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "omega",
        compressible::omegaWallFunctionFvPatchScalarField::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "k",
        compressible::kqRWallFunctionFvPatchField<scalar>::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "q",
        compressible::kqRWallFunctionFvPatchField<scalar>::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "R",
        compressible::kqRWallFunctionFvPatchField<symmTensor>::
            typeName,
        "(0 0 0 0 0 0)"
    );
}


void updateIncompressibleCase(const fvMesh& mesh)
{
    Info<< "Case treated as incompressible" << nl << endl;
    createVolScalarField(mesh, "nut", dimArea/dimTime);

    replaceBoundaryType
    (
        mesh,
        "nut",
        incompressible::nutkWallFunctionFvPatchScalarField::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "epsilon",
        incompressible::epsilonWallFunctionFvPatchScalarField::
            typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "omega",
        incompressible::omegaWallFunctionFvPatchScalarField::
            typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "k",
        incompressible::kqRWallFunctionFvPatchField<scalar>::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "q",
        incompressible::kqRWallFunctionFvPatchField<scalar>::typeName,
        "0"
    );
    replaceBoundaryType
    (
        mesh,
        "R",
        incompressible::kqRWallFunctionFvPatchField<symmTensor>::typeName,
        "(0 0 0 0 0 0)"
    );
}


int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    argList::addBoolOption
    (
        "compressible",
        "force use of compressible wall functions. Default is auto-detect."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const bool compressible = args.optionFound("compressible");

    Info<< "Updating turbulence fields to operate using new run time "
        << "selectable" << nl << "wall functions"
        << nl << endl;

    if (compressible || caseIsCompressible(mesh))
    {
        updateCompressibleCase(mesh);
    }
    else
    {
        updateIncompressibleCase(mesh);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
