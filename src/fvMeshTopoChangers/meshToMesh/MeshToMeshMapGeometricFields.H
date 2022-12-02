/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
void evaluateConstraintTypes(VolField<Type>& fld)
{
    typename VolField<Type>::
        Boundary& fldBf = fld.boundaryFieldRef();

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(Pstream::defaultCommsType);
            }
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(Pstream::defaultCommsType);
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (patchSchedule[patchEvali].init)
                {
                    tgtField.initEvaluate(Pstream::commsTypes::scheduled);
                }
                else
                {
                    tgtField.evaluate(Pstream::commsTypes::scheduled);
                }
            }
        }
    }
}


template<class Type>
void MeshToMeshMapVolFields
(
    const fvMesh& mesh,
    const fvMeshToFvMesh& mapper
)
{
    HashTable<const VolField<Type>*> fields
    (
        mesh.objectRegistry::template
            lookupClass<VolField<Type>>()
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
        mesh.objectRegistry::template
            lookupClass<VolField<Type>>()
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
            const_cast<VolField<Type>&>
            (*fieldIter());

        if (fvMeshToFvMesh::debug)
        {
            Info<< "Mapping " << field.typeName << ' ' << field.name()
                << endl;
        }

        field.reset(mapper.mapSrcToTgt(field));
        evaluateConstraintTypes(field);
        field.instance() = field.time().name();
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class PatchMapper
>
void NaNGeometricFields
(
    const fvMesh& mesh,
    const fvMeshToFvMesh& mapper
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GField;

    //- Mapper for sizing only - does not do any actual mapping.
    class patchFieldResizeMapper
    :
        public PatchMapper,
        public setSizeFieldMapper
    {
    public:

        // Constructors

            //- Construct given size
            patchFieldResizeMapper(const label size)
            :
                setSizeFieldMapper(size)
            {}
    };

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
                pf.autoMap(patchFieldResizeMapper(pf.patch().size()));
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
