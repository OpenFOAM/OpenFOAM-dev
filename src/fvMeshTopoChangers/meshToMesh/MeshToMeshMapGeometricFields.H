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
#include "meshToMesh.H"
#include "fvPatchFieldMapper.H"
#include "pointPatchFieldMapper.H"
#include "setSizeFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void evaluateConstraintTypes(GeometricField<Type, fvPatchField, volMesh>& fld)
{
    typename GeometricField<Type, fvPatchField, volMesh>::
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
    const meshToMesh& mapper
)
{
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mapper.srcRegion().objectRegistry::template
            lookupClass<GeometricField<Type, fvPatchField, volMesh>>()
    );

    // Deleted old time fields
    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        field.clearOldTimes();
    }

    fields =
    (
        mapper.srcRegion().objectRegistry::template
            lookupClass<GeometricField<Type, fvPatchField, volMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        if (meshToMesh::debug)
        {
            Info<< "Mapping " << field.typeName << ' ' << field.name()
                << endl;
        }

        field.reset(mapper.mapSrcToTgt(field));
        evaluateConstraintTypes(field);
        field.instance() = field.time().timeName();
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
    const meshToMesh& mapper
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> Gfield;

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

    HashTable<const Gfield*> fields
    (
        mapper.srcRegion().objectRegistry::template lookupClass<Gfield>()
    );

    // Deleted old time fields
    for
    (
        typename HashTable<const Gfield*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        Gfield& field = const_cast<Gfield&>(*fieldIter());
        field.clearOldTimes();
    }

    fields =
        mapper.srcRegion().objectRegistry::template lookupClass<Gfield>();

    Type NaN;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        setComponent(NaN, cmpt) = std::numeric_limits<scalar>::signaling_NaN();
    }

    for
    (
        typename HashTable<const Gfield*>::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        Gfield& field = const_cast<Gfield&>(*fieldIter());

        if (meshToMesh::debug)
        {
            Info<< "Setting to NaN " << field.typeName << ' ' << field.name()
                << endl;
        }

        const typename Gfield::Mesh& mesh = field.mesh();

        field.primitiveFieldRef().setSize(GeoMesh::size(mesh));
        field.primitiveFieldRef() = NaN;

        forAll(mesh.boundary(), patchi)
        {
            typename Gfield::Patch& pf = field.boundaryFieldRef()[patchi];

            pf.autoMap(patchFieldResizeMapper(pf.patch().size()));
            pf == NaN;
        }

        field.instance() = field.time().timeName();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
