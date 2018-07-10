/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "fvFieldDecomposer.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"
#include "processorCyclicFvPatchField.H"
#include "processorCyclicFvsPatchField.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const bool allowUnknownPatchFields
) const
{
    // 1. Create the complete field with dummy patch fields
    PtrList<fvPatchField<Type>> patchFields(boundaryAddressing_.size());

    forAll(boundaryAddressing_, patchi)
    {
        patchFields.set
        (
            patchi,
            fvPatchField<Type>::New
            (
                calculatedFvPatchField<Type>::typeName,
                procMesh_.boundary()[patchi],
                DimensionedField<Type, volMesh>::null()
            )
        );
    }

    // Create the field for the processor
    tmp<GeometricField<Type, fvPatchField, volMesh>> tresF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            Field<Type>(field.primitiveField(), cellAddressing_),
            patchFields
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& resF = tresF.ref();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& bf = resF.boundaryFieldRef();

    forAll(bf, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    resF(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorCyclicFvPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    resF(),
                    Field<Type>
                    (
                        field.primitiveField(),
                        *processorVolPatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorFvPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    resF(),
                    Field<Type>
                    (
                        field.primitiveField(),
                        *processorVolPatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
        else if (allowUnknownPatchFields)
        {
            bf.set
            (
                patchi,
                new emptyFvPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    resF()
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unknown type." << abort(FatalError);
        }
    }

    // Create the field for the processor
    return tresF;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    labelList mapAddr
    (
        labelList::subList
        (
            faceAddressing_,
            procMesh_.nInternalFaces()
        )
    );
    forAll(mapAddr, i)
    {
        mapAddr[i] -= 1;
    }

    // Create and map the internal field values
    Field<Type> internalField
    (
        field.primitiveField(),
        mapAddr
    );

    // Problem with addressing when a processor patch picks up both internal
    // faces and faces from cyclic boundaries. This is a bit of a hack, but
    // I cannot find a better solution without making the internal storage
    // mechanism for surfaceFields correspond to the one of faces in polyMesh
    // (i.e. using slices)
    Field<Type> allFaceField(field.mesh().nFaces());

    forAll(field.primitiveField(), i)
    {
        allFaceField[i] = field.primitiveField()[i];
    }

    forAll(field.boundaryField(), patchi)
    {
        const Field<Type> & p = field.boundaryField()[patchi];

        const label patchStart = field.mesh().boundaryMesh()[patchi].start();

        forAll(p, i)
        {
            allFaceField[patchStart + i] = p[i];
        }
    }


    // 1. Create the complete field with dummy patch fields
    PtrList<fvsPatchField<Type>> patchFields(boundaryAddressing_.size());

    forAll(boundaryAddressing_, patchi)
    {
        patchFields.set
        (
            patchi,
            fvsPatchField<Type>::New
            (
                calculatedFvsPatchField<Type>::typeName,
                procMesh_.boundary()[patchi],
                DimensionedField<Type, surfaceMesh>::null()
            )
        );
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tresF
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            Field<Type>(field.primitiveField(), mapAddr),
            patchFields
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& resF = tresF.ref();


    // 2. Change the fvsPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        Boundary& bf = resF.boundaryFieldRef();

    forAll(boundaryAddressing_, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            bf.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    resF(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorCyclicFvsPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    resF(),
                    Field<Type>
                    (
                        allFaceField,
                        *processorSurfacePatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorFvsPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    resF(),
                    Field<Type>
                    (
                        allFaceField,
                        *processorSurfacePatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unknown type." << abort(FatalError);
        }
    }

    // Create the field for the processor
    return tresF;
}


template<class GeoField>
void Foam::fvFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    forAll(fields, fieldi)
    {
        decomposeField(fields[fieldi])().write();
    }
}


// ************************************************************************* //
