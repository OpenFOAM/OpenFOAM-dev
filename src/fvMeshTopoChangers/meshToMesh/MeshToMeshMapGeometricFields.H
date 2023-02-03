/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::MeshToMeshMapGeometricFields

Description
    Generic internal field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef MeshToMeshMapGeometricFields_H
#define MeshToMeshMapGeometricFields_H

#include "polyMesh.H"
#include "fvMeshToFvMesh.H"
#include "fvPatchFieldMapper.H"
#include "pointPatchFieldMapper.H"
#include "setSizeFieldMapper.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void MeshToMeshMapVolFields
(
    const fvMesh& mesh,
    const fvMeshToFvMesh& mapper
)
{
    HashTable<const VolField<Type>*> fields
    (
        mesh.objectRegistry::template lookupClass<VolField<Type>>()
    );

    // Deleted old time fields
    for
    (
        typename HashTable<const VolField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        VolField<Type>& field =
            const_cast<VolField<Type>&>
            (*fieldIter());

        field.clearOldTimes();
    }

    fields =
    (
        mesh.objectRegistry::template lookupClass<VolField<Type>>()
    );

    for
    (
        typename HashTable<const VolField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        VolField<Type>& field =
            const_cast<VolField<Type>&>(*fieldIter());

        if (fvMeshToFvMesh::debug)
        {
            Info<< "Mapping " << field.typeName << ' ' << field.name()
                << endl;
        }

        field.reset(mapper.srcToTgt(field));

        field.instance() = field.time().name();
    }
}


template<class Type>
void MeshToMeshMapVolInternalFields
(
    const fvMesh& mesh,
    const fvMeshToFvMesh& mapper
)
{
    // Note: use strict flag on lookupClass to avoid picking up volFields
    HashTable<const VolInternalField<Type>*> fields
    (
        mesh.objectRegistry::template lookupClass<VolInternalField<Type>>(true)
    );

    for
    (
        typename HashTable<const VolInternalField<Type>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        VolInternalField<Type>& field =
            const_cast<VolInternalField<Type>&>(*fieldIter());

        if (fvMeshToFvMesh::debug)
        {
            Info<< "Mapping " << field.typeName << ' ' << field.name()
                << endl;
        }

        field.reset(mapper.srcToTgt<Type>(field));

        field.instance() = field.time().name();
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class SetSizePatchFieldMapper
>
void NaNGeometricFields
(
    const fvMesh& mesh,
    const fvMeshToFvMesh& mapper
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GField;

    HashTable<const GField*> fields
    (
        mesh.objectRegistry::template lookupClass<GField>()
    );

    // Deleted old time fields
    for
    (
        typename HashTable<const GField*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GField& field = const_cast<GField&>(*fieldIter());
        field.clearOldTimes();
    }

    fields = mesh.objectRegistry::template lookupClass<GField>();

    Type NaN;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        setComponent(NaN, cmpt) = std::numeric_limits<scalar>::signaling_NaN();
    }

    for
    (
        typename HashTable<const GField*>::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GField& field = const_cast<GField&>(*fieldIter());

        if (fvMeshToFvMesh::debug)
        {
            Info<< "Setting to NaN " << field.typeName << ' ' << field.name()
                << endl;
        }

        const typename GField::Mesh& mesh = field.mesh();

        field.primitiveFieldRef().setSize(GeoMesh::size(mesh));
        field.primitiveFieldRef() = NaN;

        field.boundaryFieldRef().setSize(mesh.boundary().size());

        forAll(mesh.boundary(), patchi)
        {
            if (isA<processorPolyPatch>(mesh().boundaryMesh()[patchi]))
            {
                field.boundaryFieldRef().set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        calculatedFvPatchField<Type>::typeName,
                        mesh.boundary()[patchi],
                        field
                    )
                );
            }
            else
            {
                typename GField::Patch& pf = field.boundaryFieldRef()[patchi];
                pf.map(pf, SetSizePatchFieldMapper(pf.patch().size()));
            }

            field.boundaryFieldRef()[patchi] == NaN;
        }

        field.instance() = field.time().name();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
