/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "readFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::readFields
(
    const typename GeoMesh::Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    const bool readOldTime
)
{
    // Search list of objects for fields
    IOobjectList fieldObjects(objects.lookupClass
    (
        GeometricField<Type, PatchField, GeoMesh>::typeName)
    );

    // Remove the cellProc field
    IOobjectList::iterator cellProcIter = fieldObjects.find("cellProc");
    if (cellProcIter != fieldObjects.end())
    {
        fieldObjects.erase(cellProcIter);
    }

    // Get sorted set of names (different processors might read objects in
    // different order)
    const wordList masterNames(fieldObjects.sortedNames());

    // Construct the fields
    fields.setSize(masterNames.size());

    forAll(masterNames, i)
    {
        const IOobject& io = *fieldObjects[masterNames[i]];

        fields.set
        (
            i,
            new GeometricField<Type, PatchField, GeoMesh>
            (
                io,
                mesh,
                readOldTime
            )
        );
    }
}


template<class Mesh, class GeoField>
void Foam::readFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields
)
{
    // Search list of objects for fields of type GeomField
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Construct the fields
    fields.setSize(fieldObjects.size());

    // Get sorted set of names (different processors might read objects in
    // different order)
    const wordList masterNames(fieldObjects.sortedNames());

    // Construct the fields
    fields.setSize(masterNames.size());

    forAll(masterNames, i)
    {
        const IOobject& io = *fieldObjects[masterNames[i]];

        fields.set(i, new GeoField(io, mesh));
    }
}


// ************************************************************************* //
