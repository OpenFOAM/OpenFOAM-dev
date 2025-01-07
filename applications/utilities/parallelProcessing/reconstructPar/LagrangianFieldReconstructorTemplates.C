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

#include "LagrangianFieldReconstructor.H"
#include "LagrangianFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
bool Foam::LagrangianFieldReconstructor::reconstructs
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    IOobjectList fields = objects.lookupClass(GeoField::typeName);

    if (fields.size() && selectedFields.empty())
    {
        return true;
    }

    forAllConstIter(IOobjectList, fields, fieldIter)
    {
        if (selectedFields.found(fieldIter()->name()))
        {
            return true;
        }
    }

    return false;
}


template
<
    class Type,
    template<class> class PrimitiveField,
    template<class, class, template<class> class> class GeoField
>
Foam::tmp<PrimitiveField<Type>>
Foam::LagrangianFieldReconstructor::reconstructLagrangianPrimitiveField
(
    const PtrList<GeoField<Type, LagrangianMesh, PrimitiveField>>&
        procFields
) const
{
    tmp<PrimitiveField<Type>> tresult
    (
        new PrimitiveField<Type>(completeMesh_.size())
    );
    PrimitiveField<Type>& result = tresult.ref();

    // Expand dynamic primitive fields to their full size
    result.setSize(completeMesh_.size());

    label i0 = 0;
    forAll(procMeshes_, proci)
    {
        SubList<Type>(result, procMeshes_[proci].size(), i0) =
            procFields[proci].primitiveField();

        i0 += procMeshes_[proci].size();
    }

    return tresult;
}


template<class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, Foam::LagrangianMesh, PrimitiveField>>
Foam::LagrangianFieldReconstructor::reconstructLagrangianField
(
    const PtrList<DimensionedField<Type, LagrangianMesh, PrimitiveField>>&
        procFields
) const
{
    return
        tmp<DimensionedField<Type, LagrangianMesh, PrimitiveField>>
        (
            new DimensionedField<Type, LagrangianMesh, PrimitiveField>
            (
                IOobject
                (
                    procFields[0].name(),
                    completeMesh_.time().name(),
                    completeMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                completeMesh_,
                procFields[0].dimensions(),
                reconstructLagrangianPrimitiveField(procFields)()
            )
        );
}


template<class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::GeometricField<Type, Foam::LagrangianMesh, PrimitiveField>>
Foam::LagrangianFieldReconstructor::reconstructLagrangianField
(
    const PtrList<GeometricField<Type, LagrangianMesh, PrimitiveField>>&
        procFields
) const
{
    // Construct a list of patch fields. Clone the global fields. Ignore
    // processor-specific ones. Note that we supply the processor-zero field as
    // the internal field reference. This is a placeholder and will get
    // replaced by the complete field when this set of patch fields is cloned
    // to produce the final field.
    PtrList<LagrangianPatchField<Type>> completePatchFields
    (
        completeMesh_.boundary().size()
    );
    forAll(completeMesh_.boundary(), completePatchi)
    {
        completePatchFields.set
        (
            completePatchi,
            procFields[0].boundaryField()[completePatchi].clone
            (
                procFields[0]
            )
        );
    }

    // Construct and return the complete field
    return
        tmp<GeometricField<Type, LagrangianMesh, PrimitiveField>>
        (
            new GeometricField<Type, LagrangianMesh, PrimitiveField>
            (
                IOobject
                (
                    procFields[0].name(),
                    completeMesh_.time().name(),
                    completeMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                completeMesh_,
                procFields[0].dimensions(),
                reconstructLagrangianPrimitiveField(procFields)(),
                completePatchFields,
                procFields[0].sources().table()
            )
        );
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::LagrangianFieldReconstructor::reconstructLagrangianField
(
    const IOobject& fieldIo
) const
{
    // Read the fields from all the processors
    PtrList<GeoField> procFields(procMeshes_.size());
    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeoField
            (
                IOobject
                (
                    fieldIo.name(),
                    procMeshes_[proci].time().name(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci]
            )
        );
    }

    // Reconstruct and return
    return reconstructLagrangianField(procFields);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
void Foam::LagrangianFieldReconstructor::reconstructFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    IOobjectList fields = objects.lookupClass(GeoField::typeName);

    if (fields.size())
    {
        Info<< nl << "    Reconstructing " << GeoField::typeName << "s"
            << nl << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructLagrangianField<GeoField>(*fieldIter())().write();
            }
        }
    }
}


// ************************************************************************* //
