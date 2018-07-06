/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "surfaceInterpolate.H"
#include "volFields.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::surfaceInterpolate::interpolateFields
(
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    // Convert field to map
    HashTable<word> fieldMap(2*fieldSet_.size());
    forAll(fieldSet_, i)
    {
        fieldMap.insert(fieldSet_[i].first(), fieldSet_[i].second());
    }


    HashTable<const VolFieldType*> flds(obr_.lookupClass<VolFieldType>());

    forAllConstIter(typename HashTable<const VolFieldType*>, flds, iter)
    {
        const VolFieldType& fld = *iter();

        if (fieldMap.found(fld.name()))
        {
            // const word sName = "interpolate(" + fld.name() + ')';
            const word& sName = fieldMap[fld.name()];

            if (obr_.found(sName))
            {
                Info<< "        surface field " << sName << " already exists"
                    << endl;
            }
            else
            {
                label sz = sflds.size();
                sflds.setSize(sz+1);
                sflds.set
                (
                    sz,
                    new SurfaceFieldType(sName, linearInterpolate(fld))
                );

                Info<< "        interpolated " << fld.name() << " to create "
                    << sflds[sz].name() << endl;
            }
        }
    }
}


// ************************************************************************* //
