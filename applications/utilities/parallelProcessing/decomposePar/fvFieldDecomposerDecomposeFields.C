/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
Foam::tmp<Foam::Field<Type>> Foam::fvFieldDecomposer::mapField
(
    const Field<Type>& field,
    const labelUList& mapAndSign,
    const bool applyFlip
)
{
    tmp<Field<Type>> tfld(new Field<Type>(mapAndSign.size()));
    Field<Type>& fld = tfld.ref();

    if (applyFlip)
    {
        forAll(mapAndSign, i)
        {
            if (mapAndSign[i] < 0)
            {
                fld[i] = -field[-mapAndSign[i] - 1];
            }
            else
            {
                fld[i] = field[mapAndSign[i] - 1];
            }
        }
    }
    else
    {
        // Ignore face flipping
        fld.map(field, mag(mapAndSign) - 1);
    }
    return tfld;
}


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
                IOobject::NO_WRITE,
                false
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
        const fvPatch& procPatch = procMesh_.boundary()[patchi];

        label fromPatchi = boundaryAddressing_[patchi];
        if (fromPatchi < 0 && isA<processorCyclicFvPatch>(procPatch))
        {
            const label referPatchi =
                refCast<const processorCyclicPolyPatch>
                (procPatch.patch()).referPatchID();
            if (field.boundaryField()[referPatchi].overridesConstraint())
            {
                fromPatchi = boundaryAddressing_[referPatchi];
            }
        }

        if (fromPatchi >= 0)
        {
            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    field.boundaryField()[fromPatchi],
                    procPatch,
                    resF(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procPatch))
        {
            bf.set
            (
                patchi,
                new processorCyclicFvPatchField<Type>
                (
                    procPatch,
                    resF(),
                    (*processorVolPatchFieldDecomposerPtrs_[patchi])
                    (
                        field.primitiveField()
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procPatch))
        {
            bf.set
            (
                patchi,
                new processorFvPatchField<Type>
                (
                    procPatch,
                    resF(),
                    (*processorVolPatchFieldDecomposerPtrs_[patchi])
                    (
                        field.primitiveField()
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
                    procPatch,
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
    // Apply flipping to surfaceScalarFields only
    const bool doFlip = (pTraits<Type>::nComponents == 1);


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
                IOobject::NO_WRITE,
                false
            ),
            procMesh_,
            field.dimensions(),
            mapField
            (
                field,
                labelList::subList
                (
                    faceAddressing_,
                    procMesh_.nInternalFaces()
                ),
                doFlip
            ),
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
        const fvPatch& procPatch = procMesh_.boundary()[patchi];

        if (boundaryAddressing_[patchi] >= 0)
        {
            bf.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procPatch,
                    resF(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procPatch))
        {
            // Do our own mapping. Avoids a lot of mapping complexity.
            bf.set
            (
                patchi,
                new processorCyclicFvsPatchField<Type>
                (
                    procPatch,
                    resF(),
                    mapField
                    (
                        allFaceField,
                        procPatch.patchSlice(faceAddressing_),
                        doFlip
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procPatch))
        {
            // Do our own mapping. Avoids a lot of mapping complexity.
            bf.set
            (
                patchi,
                new processorFvsPatchField<Type>
                (
                    procPatch,
                    resF(),
                    mapField
                    (
                        allFaceField,
                        procPatch.patchSlice(faceAddressing_),
                        doFlip
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
