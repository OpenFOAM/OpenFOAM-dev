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

\*---------------------------------------------------------------------------*/

#include "backwardsCompatibilityWallFunctions.H"
#include "Time.H"
#include "OSspecific.H"

#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class PatchType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
autoCreateWallFunctionField
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    IOobject nutHeader
    (
        "nut",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (nutHeader.headerOk())
    {
        return tmp<fieldType>
        (
            new fieldType
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "--> Upgrading " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        // Read existing field
        IOobject ioObj
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        tmp<fieldType> fieldOrig
        (
            new fieldType
            (
                ioObj,
                mesh
            )
        );

        // rename file
        Info<< "    Backup original " << fieldName << " to "
            << fieldName << ".old" << endl;
        mvBak(ioObj.objectPath(), "old");


        PtrList<fvPatchField<Type> > newPatchFields(mesh.boundary().size());

        forAll(newPatchFields, patchI)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchI]))
            {
                newPatchFields.set
                (
                    patchI,
                    new PatchType
                    (
                        mesh.boundary()[patchI],
                        fieldOrig().dimensionedInternalField()
                    )
                );
                newPatchFields[patchI] == fieldOrig().boundaryField()[patchI];
            }
            else
            {
                newPatchFields.set
                (
                    patchI,
                    fieldOrig().boundaryField()[patchI].clone()
                );
            }
        }

        tmp<fieldType> fieldNew
        (
            new fieldType
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
                fieldOrig().dimensions(),
                fieldOrig().internalField(),
                newPatchFields
            )
        );

        Info<< "    Writing updated " << fieldName << endl;
        fieldNew().write();

        return fieldNew;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
