/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "GeometricField.H"
#include "fvPatchField.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::helpTypes::helpBoundary::fieldConditions
(
    const IOobject& io,
    const bool write
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (io.headerClassName() == fieldType::typeName)
    {
        wordList types
        (
            fvPatchField<Type>::dictionaryConstructorTablePtr_->sortedToc()
        );

        if (write)
        {
            Info<< "Available boundary conditions for "
                << pTraits<Type>::typeName << " field: " << io.name() << nl;

            forAll(types, i)
            {
                Info<< "    " << types[i] << nl;
            }

            Info<< endl;
        }

        return types;
    }

    return wordList();
}


template<class Type>
void Foam::helpTypes::helpBoundary::fixedValueFieldConditions
(
    const IOobject& io
) const
{
    wordList types(fieldConditions<Type>(io, false));

    if (!types.size())
    {
        return;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& mesh = dynamic_cast<const fvMesh&>(io.db());

    fieldType fld
    (
        IOobject
        (
            "dummy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,   
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
    );


    Info<< "Fixed value boundary conditions for "
        << pTraits<Type>::typeName << " field: " << io.name() << nl;

    // throw exceptions to avoid fatal errors when casting from generic patch
    // type to incompatible patch type
    FatalIOError.throwExceptions();
    FatalError.throwExceptions();

    bool foundFixed = false;
    forAll(types, i)
    {
        const word& patchType = types[i];

        try
        {
            polyPatch pp
            (
                "defaultFaces",
                0,
                mesh.nInternalFaces(),
                0,
                mesh.boundaryMesh(),
                patchType
            );

            fvPatch fvp(pp, mesh.boundary());

            tmp<fvPatchField<Type> > pf
            (
                fvPatchField<Type>::New
                (
                    patchType,
                    fvp,
                    fld
                )
            );

            if (pf().fixesValue())
            {
                Info<< "    " << patchType << nl;
                foundFixed = true;
            }
        }
        catch (...)
        {
            // do nothing
        }
    }

    if (!foundFixed)
    {
        // no conditions???
        Info<< "    none" << nl;
    }

    Info<< endl;
}


// ************************************************************************* //
