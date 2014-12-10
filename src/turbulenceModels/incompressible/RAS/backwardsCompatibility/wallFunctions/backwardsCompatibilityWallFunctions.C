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

\*---------------------------------------------------------------------------*/

#include "backwardsCompatibilityWallFunctions.H"

#include "volFields.H"
#include "calculatedFvPatchField.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "nutLowReWallFunctionFvPatchScalarField.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "kqRWallFunctionFvPatchField.H"
#include "omegaWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> autoCreateNut
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    IOobject nutHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (nutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(nutHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList nutBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (isA<wallFvPatch>(bm[patchI]))
            {
                nutBoundaryTypes[patchI] =
                    nutkWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                nutBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> nut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimArea/dimTime, 0.0),
                nutBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        nut().write();

        return nut;
    }
}


tmp<volScalarField> autoCreateLowReNut
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    IOobject nutHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (nutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(nutHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList nutBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (isA<wallFvPatch>(bm[patchI]))
            {
                nutBoundaryTypes[patchI] =
                    nutLowReWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                nutBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> nut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimArea/dimTime, 0.0),
                nutBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        nut().write();

        return nut;
    }
}


tmp<volScalarField> autoCreateEpsilon
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            epsilonWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh
        );
}


tmp<volScalarField> autoCreateOmega
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            omegaWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh
        );
}


tmp<volScalarField> autoCreateK
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh
        );
}


tmp<volScalarField> autoCreateQ
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh
        );
}


tmp<volSymmTensorField> autoCreateR
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            symmTensor,
            kqRWallFunctionFvPatchField<symmTensor>
        >
        (
            fieldName,
            mesh
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

