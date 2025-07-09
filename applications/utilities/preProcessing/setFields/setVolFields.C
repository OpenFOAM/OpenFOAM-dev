/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "volFields.H"
#include "processorFvPatch.H"
#include "CompactListList.H"
#include "PtrDictionary.H"
#include "wordReListMatcher.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void setVolField
(
    const word& fieldName,
    const IOobject& fieldHeader,
    const fvMesh& mesh,
    const labelList& selectedCells,
    const PtrDictionary<wordReList>& extrapolatePatches,
    Istream& fieldValueStream
)
{
    if (fieldHeader.headerClassName() == VolField<Type>::typeName)
    {
        Info<< "    Setting "
            << VolField<Type>::typeName << "::Internal " << fieldName << endl;

        VolField<Type> field(fieldHeader, mesh);

        const Type value = pTraits<Type>(fieldValueStream);

        if (&selectedCells == &labelList::null())
        {
            field.primitiveFieldRef() = value;
        }
        else
        {
            forAll(selectedCells, celli)
            {
                field[selectedCells[celli]] = value;
            }
        }

        typename VolField<Type>::
            Boundary& fieldBf = field.boundaryFieldRef();

        forAll(field.boundaryField(), patchi)
        {
            if
            (
                extrapolatePatches.found(mesh.boundary()[patchi].name())
             && wordReListMatcher
                (
                    extrapolatePatches[mesh.boundary()[patchi].name()]
                ).match(fieldName)
            )
            {
                fieldBf[patchi] == fieldBf[patchi].patchInternalField();
            }
            else
            {
                fieldBf[patchi] = fieldBf[patchi].patchInternalField();
            }
        }

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << fieldName << endl;
        }
    }
}


template<class Type>
void setPatchField
(
    const word& fieldName,
    const IOobject& fieldHeader,
    const fvMesh& mesh,
    const labelList& selectedFaces,
    Istream& fieldValueStream
)
{
    if (fieldHeader.headerClassName() == VolField<Type>::typeName)
    {
        Info<< "    Setting "
            << VolField<Type>::typeName << "::Boundary " << fieldName << endl;

        // Read the field
        VolField<Type> field(fieldHeader, mesh);
        typename VolField<Type>::Boundary& fieldBf = field.boundaryFieldRef();

        // Read the value
        const Type value = pTraits<Type>(fieldValueStream);

        // Determine the number of non-processor patches
        label nNonProcPatches = 0;
        forAll(fieldBf, patchi)
        {
            if (isA<processorFvPatch>(mesh.boundary()[patchi]))
            {
                break;
            }
            nNonProcPatches = patchi + 1;
        }

        // Create a copy of the boundary field
        typename VolField<Type>::Boundary fieldBfCopy
        (
            VolField<Type>::Internal::null(),
            fieldBf
        );

        // Loop selected faces and set values in the copied boundary field
        bool haveWarnedInternal = false, haveWarnedProc = false;
        labelList nonProcPatchNChangedFaces(nNonProcPatches, 0);
        forAll(selectedFaces, i)
        {
            const label facei = selectedFaces[i];

            if (mesh.isInternalFace(facei))
            {
                if (!haveWarnedInternal)
                {
                    WarningInFunction
                        << "Ignoring internal face " << facei
                        << ". Suppressing further warnings." << endl;
                    haveWarnedInternal = true;
                }
            }
            else
            {
                const labelUList patches =
                    mesh.polyBFacePatches()[facei - mesh.nInternalFaces()];
                const labelUList patchFaces =
                    mesh.polyBFacePatchFaces()[facei - mesh.nInternalFaces()];

                forAll(patches, i)
                {
                    if (patches[i] >= nNonProcPatches)
                    {
                        if (!haveWarnedProc)
                        {
                            WarningInFunction
                                << "Ignoring face " << patchFaces[i]
                                << " of processor patch " << patches[i]
                                << ". Suppressing further warnings." << endl;
                            haveWarnedProc = true;
                        }
                    }
                    else
                    {
                        fieldBfCopy[patches[i]][patchFaces[i]] = value;
                        nonProcPatchNChangedFaces[patches[i]] ++;
                    }
                }
            }
        }
        Pstream::listCombineGather
        (
            nonProcPatchNChangedFaces,
            plusEqOp<label>()
        );
        Pstream::listCombineScatter
        (
            nonProcPatchNChangedFaces
        );

        // Reassign boundary values
        forAll(nonProcPatchNChangedFaces, patchi)
        {
            if (nonProcPatchNChangedFaces[patchi] > 0)
            {
                Info<< "    On patch "
                    << field.boundaryField()[patchi].patch().name()
                    << " set " << nonProcPatchNChangedFaces[patchi]
                    << " values" << endl;
                fieldBf[patchi] == fieldBfCopy[patchi];
            }
        }

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << exit(FatalError);
        }
    }
}


void setVolFields
(
    const fvMesh& mesh,
    const dictionary& fieldsDict,
    const labelList& selectedCells,
    const PtrDictionary<wordReList>& extrapolatePatches
)
{
    forAllConstIter(dictionary, fieldsDict, iter)
    {
        const word& fieldName = iter().keyword();

        // Check the current time directory
        IOobject fieldHeader
        (
            fieldName,
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ
        );

        // Check the "constant" directory
        if (!fieldHeader.headerOk())
        {
            fieldHeader = IOobject
            (
                fieldName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            );
        }

        // Check field exists
        if (fieldHeader.headerOk())
        {
            #define SetVolField(Type, nullArg)                                 \
                setVolField<Type>                                              \
                (                                                              \
                    fieldName,                                                 \
                    fieldHeader,                                               \
                    mesh,                                                      \
                    selectedCells,                                             \
                    extrapolatePatches,                                        \
                    iter().stream()                                            \
                );

            FOR_ALL_FIELD_TYPES(SetVolField);

            #undef SetCellFieldType
        }
        else
        {
            WarningInFunction
                << "Field " << fieldName << " not found" << endl;
        }
    }
}


void setPatchFields
(
    const fvMesh& mesh,
    const dictionary& fieldsDict,
    const labelList& selectedCells
)
{
    forAllConstIter(dictionary, fieldsDict, iter)
    {
        const word& fieldName = iter().keyword();

        // Check the current time directory
        IOobject fieldHeader
        (
            fieldName,
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ
        );

        // Check the "constant" directory
        if (!fieldHeader.headerOk())
        {
            fieldHeader = IOobject
            (
                fieldName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            );
        }

        // Check field exists
        if (fieldHeader.headerOk())
        {
            #define SetPatchField(Type, nullArg)                               \
                setPatchField<Type>                                            \
                (                                                              \
                    fieldName,                                                 \
                    fieldHeader,                                               \
                    mesh,                                                      \
                    selectedCells,                                             \
                    iter().stream()                                            \
                );

            FOR_ALL_FIELD_TYPES(SetPatchField);

            #undef SetCellFieldType
        }
        else
        {
            WarningInFunction
                << "Field " << fieldName << " not found" << endl;
        }
    }
}


// ************************************************************************* //
