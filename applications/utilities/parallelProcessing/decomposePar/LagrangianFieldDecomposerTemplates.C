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

#include "LagrangianFieldDecomposer.H"
#include "LagrangianFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
bool Foam::LagrangianFieldDecomposer::decomposes(const IOobjectList& objects)
{
    return !objects.lookupClass(GeoField::typeName).empty();
}


template<class Type, template<class> class PrimitiveField>
Foam::PtrList
<
    Foam::DimensionedField<Type, Foam::LagrangianMesh, PrimitiveField>
>
Foam::LagrangianFieldDecomposer::decomposeLagrangianField
(
    const DimensionedField<Type, LagrangianMesh, PrimitiveField>& completeField
) const
{
    PtrList<DimensionedField<Type, LagrangianMesh, PrimitiveField>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new DimensionedField<Type, LagrangianMesh, PrimitiveField>
            (
                IOobject
                (
                    completeField.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci],
                completeField.dimensions(),
                PrimitiveField<Type>
                (
                    completeField,
                    particleProcAddressing_[proci]
                )
            )
        );
    }

    return procFields;
}


template<class Type, template<class> class PrimitiveField>
Foam::PtrList
<
    Foam::GeometricField<Type, Foam::LagrangianMesh, PrimitiveField>
>
Foam::LagrangianFieldDecomposer::decomposeLagrangianField
(
    const GeometricField<Type, LagrangianMesh, PrimitiveField>& completeField
) const
{
    PtrList<GeometricField<Type, LagrangianMesh, PrimitiveField>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        // Construct a list of patch fields. Clone the global fields.
        // Null-construct new processor-specific ones. Note that we supply the
        // complete field as the internal field reference. This is a
        // placeholder and will get replaced with the processor field when this
        // set of patch fields is cloned to produce the final field.
        PtrList<LagrangianPatchField<Type>> procPatchFields
        (
            procMeshes_[proci].boundary().size()
        );
        forAll(completeField.boundaryField(), completePatchi)
        {
            procPatchFields.set
            (
                completePatchi,
                completeField.boundaryField()[completePatchi].clone
                (
                    completeField
                )
            );
        }
        for
        (
            label procPatchi = completeMesh_.boundary().size();
            procPatchi < procMeshes_[proci].boundary().size();
            ++ procPatchi
        )
        {
            procPatchFields.set
            (
                procPatchi,
                LagrangianPatchField<Type>::New
                (
                    procMeshes_[proci].boundary()[procPatchi].type(),
                    procMeshes_[proci].boundary()[procPatchi],
                    completeField
                )
            );
        }

        // Construct the field
        procFields.set
        (
            proci,
            new GeometricField<Type, LagrangianMesh, PrimitiveField>
            (
                IOobject
                (
                    completeField.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci],
                completeField.dimensions(),
                PrimitiveField<Type>
                (
                    completeField,
                    particleProcAddressing_[proci]
                ),
                procPatchFields,
                completeField.sources().table()
            )
        );
    }

    return procFields;
}


template<class GeoField>
Foam::PtrList<GeoField>
Foam::LagrangianFieldDecomposer::decomposeLagrangianField
(
    const IOobject& fieldIo
) const
{
    return
        decomposeLagrangianField
        (
            GeoField
            (
                IOobject
                (
                    fieldIo.name(),
                    completeMesh_.time().name(),
                    completeMesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                completeMesh_
            )
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
void Foam::LagrangianFieldDecomposer::decomposeFields
(
    const IOobjectList& objects
) const
{
    IOobjectList fields = objects.lookupClass(GeoField::typeName);

    if (fields.size())
    {
        Info<< nl << "    Decomposing " << GeoField::typeName << "s"
            << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            Info<< "        " << fieldIter()->name() << endl;

            PtrList<GeoField> procFields =
                decomposeLagrangianField<GeoField>(*fieldIter());

            forAll(procFields, proci)
            {
                procFields[proci].write();
            }
        }
    }
}


// ************************************************************************* //
