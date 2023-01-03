/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "SlicedGeometricField.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::tmp<Foam::FieldField<GeoMesh::template PatchField, Type>>
Foam::SlicedGeometricField<Type, GeoMesh>::slicedBoundaryField
(
    const Mesh& mesh,
    const Field<Type>& completeField,
    const bool preserveCouples,
    const bool preserveProcessorOnly
)
{
    tmp<FieldField<GeoMesh::template PatchField, Type>> tbf
    (
        new FieldField<GeoMesh::template PatchField, Type>
        (
            mesh.boundary().size()
        )
    );
    FieldField<GeoMesh::template PatchField, Type>& bf = tbf.ref();

    forAll(mesh.boundary(), patchi)
    {
        if
        (
            preserveCouples
         && mesh.boundary()[patchi].coupled()
         && (
               !preserveProcessorOnly
            || isA<processorFvPatch>(mesh.boundary()[patchi])
            )
        )
        {
            // For coupled patched construct the correct patch field type
            bf.set
            (
                patchi,
                Patch::New
                (
                    mesh.boundary()[patchi].type(),
                    mesh.boundary()[patchi],
                    *this
                )
            );

            // Initialise the values on the coupled patch to those of the slice
            // of the given field.
            // Note: these will usually be over-ridden by the boundary field
            // evaluation e.g. in the case of processor and cyclic patches.
            bf[patchi] = SlicedPatch
            (
                mesh.boundary()[patchi],
                DimensionedField<Type, GeoMesh>::null(),
                completeField
            );
        }
        else
        {
            bf.set
            (
                patchi,
                new SlicedPatch
                (
                    mesh.boundary()[patchi],
                    DimensionedField<Type, GeoMesh>::null(),
                    completeField
                )
            );
        }
    }

    return tbf;
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::FieldField<GeoMesh::template PatchField, Type>>
Foam::SlicedGeometricField<Type, GeoMesh>::slicedBoundaryField
(
    const Mesh& mesh,
    const FieldField<GeoMesh::template PatchField, Type>& bField,
    const bool preserveCouples
)
{
    tmp<FieldField<GeoMesh::template PatchField, Type>> tbf
    (
        new FieldField<GeoMesh::template PatchField, Type>
        (
            mesh.boundary().size()
        )
    );
    FieldField<GeoMesh::template PatchField, Type>& bf = tbf.ref();

    forAll(mesh.boundary(), patchi)
    {
        if (preserveCouples && mesh.boundary()[patchi].coupled())
        {
            // For coupled patched construct the correct patch field type
            bf.set
            (
                patchi,
                Patch::New
                (
                    mesh.boundary()[patchi].type(),
                    mesh.boundary()[patchi],
                    *this
                )
            );

            // Assign field
            bf[patchi] == bField[patchi];
        }
        else
        {
            // Create unallocated copy of patch field
            bf.set
            (
                patchi,
                new SlicedPatch
                (
                    mesh.boundary()[patchi],
                    DimensionedField<Type, GeoMesh>::null(),
                    bField[patchi]
                )
            );
        }
    }

    return tbf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::SlicedGeometricField<Type, GeoMesh>::SlicedGeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& completeField,
    const bool preserveCouples
)
:
    GeometricField<Type, GeoMesh>
    (
        io,
        mesh,
        ds,
        Field<Type>(),
        slicedBoundaryField(mesh, completeField, preserveCouples)
    )
{
    // Set the internalField to the slice of the complete field
    UList<Type>::shallowCopy
    (
        typename Field<Type>::subField(completeField, GeoMesh::size(mesh))
    );

    correctBoundaryConditions();
}


template<class Type, class GeoMesh>
Foam::SlicedGeometricField<Type, GeoMesh>::SlicedGeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& completeIField,
    const Field<Type>& completeBField,
    const bool preserveCouples,
    const bool preserveProcessorOnly
)
:
    GeometricField<Type, GeoMesh>
    (
        io,
        mesh,
        ds,
        Field<Type>(),
        slicedBoundaryField
        (
            mesh,
            completeBField,
            preserveCouples,
            preserveProcessorOnly
        )
    )
{
    // Set the internalField to the slice of the complete field
    UList<Type>::shallowCopy
    (
        typename Field<Type>::subField(completeIField, GeoMesh::size(mesh))
    );

    correctBoundaryConditions();
}


template<class Type, class GeoMesh>
Foam::SlicedGeometricField<Type, GeoMesh>::SlicedGeometricField
(
    const IOobject& io,
    const GeometricField<Type, GeoMesh>& gf,
    const bool preserveCouples
)
:
    GeometricField<Type, GeoMesh>
    (
        io,
        gf.mesh(),
        gf.dimensions(),
        Field<Type>(),
        slicedBoundaryField(gf.mesh(), gf.boundaryField(), preserveCouples)
    )
{
    // Set the internalField to the supplied internal field
    UList<Type>::shallowCopy(gf.primitiveField());

    correctBoundaryConditions();
}


template<class Type, class GeoMesh>
Foam::SlicedGeometricField<Type, GeoMesh>::SlicedGeometricField
(
    const SlicedGeometricField<Type, GeoMesh>& gf
)
:
    GeometricField<Type, GeoMesh>
    (
        gf,
        gf.mesh(),
        gf.dimensions(),
        Field<Type>(),
        slicedBoundaryField(gf.mesh(), gf.boundaryField(), true)
    )
{
    // Set the internalField to the supplied internal field
    UList<Type>::shallowCopy(gf.primitiveField());
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::SlicedGeometricField<Type, GeoMesh>>
Foam::SlicedGeometricField<Type, GeoMesh>::clone() const
{
    return tmp<SlicedGeometricField<Type, GeoMesh>>
    (
        new SlicedGeometricField<Type, GeoMesh>(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::SlicedGeometricField<Type, GeoMesh>::~SlicedGeometricField()
{
    // Set the internalField storage pointer to nullptr before its destruction
    // to protect the field it a slice of.
    UList<Type>::shallowCopy(UList<Type>(nullptr, 0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::tmp<Foam::Field<Type>>
Foam::SlicedGeometricField<Type, GeoMesh>::splice() const
{
    const Mesh& mesh = this->mesh();

    label completeSize = GeoMesh::size(mesh);

    forAll(mesh.boundaryMesh(), patchi)
    {
        completeSize += mesh.boundaryMesh()[patchi].size();
    }

    tmp<Field<Type>> tCompleteField(new Field<Type>(completeSize));
    Field<Type>& completeField(tCompleteField.ref());

    typename Field<Type>::subField(completeField, GeoMesh::size(mesh))
        = this->primitiveField();

    label start = GeoMesh::size(mesh);

    forAll(mesh.boundaryMesh(), patchi)
    {
        if
        (
            mesh.boundary()[patchi].size()
         == mesh.boundaryMesh()[patchi].size()
        )
        {
            typename Field<Type>::subField
            (
                completeField,
                mesh.boundary()[patchi].size(),
                start
            ) = this->boundaryField()[patchi];
        }
        else
        {
            typename Field<Type>::subField
            (
                completeField,
                mesh.boundaryMesh()[patchi].size(),
                start
            ) = Zero;
        }

        start += mesh.boundaryMesh()[patchi].size();
    }

    return tCompleteField;
}


template<class Type, class GeoMesh>
void Foam::SlicedGeometricField<Type, GeoMesh>::correctBoundaryConditions()
{
    GeometricField<Type, GeoMesh>::correctBoundaryConditions();
}


// ************************************************************************* //
