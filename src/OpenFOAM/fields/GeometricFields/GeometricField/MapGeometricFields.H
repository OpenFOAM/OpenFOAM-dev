/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::MapInternalField

Description
    Generic internal field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef MapGeometricFields_H
#define MapGeometricFields_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper, class GeoMesh>
class MapInternalField
{
public:

    MapInternalField()
    {}

    void operator()
    (
        Field<Type>& field,
        const MeshMapper& mapper
    ) const;
};


//- Generic Geometric field mapper.
//  For "real" mapping, add template specialisations
//  for mapping of internal fields depending on mesh type.
template
<
    class Type,
    template<class> class PatchField,
    class MeshMapper,
    class GeoMesh
>
void MapGeometricFields
(
    const MeshMapper& mapper
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mapper.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        // Note: check can be removed once pointFields are actually stored on
        //      the pointMesh instead of now on the polyMesh!
        if (&field.mesh() == &mapper.mesh())
        {
            field.storeOldTimes();
        }
    }

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        if (&field.mesh() == &mapper.mesh())
        {
            if (polyMesh::debug)
            {
                Info<< "Mapping " << field.typeName << ' ' << field.name()
                    << endl;
            }

            // Map the internal field
            MapInternalField<Type, MeshMapper, GeoMesh>()
            (
                field.primitiveFieldRef(),
                mapper
            );

            // Map the patch fields
            typename GeometricField<Type, PatchField, GeoMesh>
            ::Boundary& bfield = field.boundaryFieldRef();
            forAll(bfield, patchi)
            {
                // Cannot check sizes for patch fields because of
                // empty fields in FV and because point fields get their size
                // from the patch which has already been resized
                //

                bfield[patchi].autoMap(mapper.boundaryMap()[patchi]);
            }

            field.instance() = field.time().timeName();
        }
        else if (polyMesh::debug)
        {
            Info<< "Not mapping " << field.typeName << ' ' << field.name()
                << " since originating mesh differs from that of mapper."
                << endl;
        }
    }
}


template<class GeoField>
void AddPatchFields
(
    objectRegistry& obr,
    const label insertPatchi,
    const dictionary& patchFieldDict,
    const word& defaultPatchFieldType,
    const typename GeoField::value_type& defaultPatchValue
)
{
    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());

    const wordList fldNames(flds.sortedToc());

    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        if (bfld.size() != fld.mesh().boundary().size())
        {
            FatalErrorInFunction << "bfld size:" << bfld.size()
                << " mesh size:" << fld.mesh().boundary().size()
                << exit(FatalError);
        }

        if (patchFieldDict.found(fld.name()))
        {
            bfld.set
            (
                insertPatchi,
                GeoField::Patch::New
                (
                    fld.mesh().boundary()[insertPatchi],
                    fld(),
                    patchFieldDict.subDict(fld.name())
                )
            );
        }
        else
        {
            bfld.set
            (
                insertPatchi,
                GeoField::Patch::New
                (
                    defaultPatchFieldType,
                    fld.mesh().boundary()[insertPatchi],
                    fld()
                )
            );
            bfld[insertPatchi] == defaultPatchValue;
        }
    }
}


template<class GeoField>
void ReorderPatchFields
(
    objectRegistry& obr,
    const labelUList& newToOld
)
{
    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());
    const wordList fldNames(flds.sortedToc());
    forAll(fldNames, i)
    {
        GeoField& fld = *flds[fldNames[i]];
        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        bfld.shuffle(newToOld);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
