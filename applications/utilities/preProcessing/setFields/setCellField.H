/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

Class
    setCellField

Description

\*---------------------------------------------------------------------------*/

#ifndef setCellField_H
#define setCellField_H

#include "processorFvPatch.H"
#include "CompactListList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool setCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream
)
{
    if (fieldTypeDesc != VolField<Type>::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    typeIOobject<VolField<Type>> fieldHeader
    (
        fieldName,
        mesh.time().name(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.headerOk())
    {
        fieldHeader = typeIOobject<VolField<Type>>
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
        Info<< "    Setting internal values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

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
            fieldBf[patchi] = fieldBf[patchi].patchInternalField();
        }

        if (!field.write())
        {
            FatalErrorInFunction
              << "Failed writing field " << fieldName << endl;
        }
    }
    else
    {
        WarningInFunction
          << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setCellField
{

public:

    setCellField()
    {}

    autoPtr<setCellField> clone() const
    {
        return autoPtr<setCellField>(new setCellField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(selectedCells)
        {}

        autoPtr<setCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setCellFieldType<scalar>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<vector>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<symmTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<tensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setCellField>(new setCellField());
        }
    };
};


template<class Type>
bool setFaceFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedFaces,
    Istream& fieldValueStream
)
{
    if
    (
        fieldTypeDesc != word::null
     && fieldTypeDesc != VolField<Type>::typeName + "Value"
    )
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    typeIOobject<VolField<Type>> fieldHeader
    (
        fieldName,
        mesh.time().name(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.headerOk())
    {
        fieldHeader = typeIOobject<VolField<Type>>
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
        Info<< "    Setting patchField values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

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
    else
    {
        WarningInFunction
          << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setFaceField
{

public:

    setFaceField()
    {}

    autoPtr<setFaceField> clone() const
    {
        return autoPtr<setFaceField>(new setFaceField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedFaces_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedFaces)
        :
            mesh_(mesh),
            selectedFaces_(selectedFaces)
        {}

        autoPtr<setFaceField> operator()(Istream& fieldStream) const
        {
            word fieldType(fieldStream);

            if
            (
               !(
                    setFaceFieldType<scalar>
                        (fieldType, mesh_, selectedFaces_, fieldStream)
                 || setFaceFieldType<vector>
                        (fieldType, mesh_, selectedFaces_, fieldStream)
                 || setFaceFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedFaces_, fieldStream)
                 || setFaceFieldType<symmTensor>
                        (fieldType, mesh_, selectedFaces_, fieldStream)
                 || setFaceFieldType<tensor>
                        (fieldType, mesh_, selectedFaces_, fieldStream)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setFaceField>(new setFaceField());
        }
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
