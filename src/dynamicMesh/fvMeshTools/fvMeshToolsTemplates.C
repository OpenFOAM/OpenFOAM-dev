/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "fvMeshTools.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    typename GeoField::Mesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    objectRegistry& obr = const_cast<objectRegistry&>(mesh.thisDb());

    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        if
        (
            patchFieldDict.found(fld.name())
        || !fvPatch::constraintType(mesh.boundary()[patchi].type())
        )
        {
            bfld.set
            (
                patchi,
                GeoField::Patch::New
                (
                    mesh.boundary()[patchi],
                    fld(),
                    patchFieldDict.subDict(fld.name())
                )
            );
        }
    }
}


template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    typename GeoField::Mesh& mesh,
    const label patchi,
    const typename GeoField::value_type& value
)
{
    objectRegistry& obr = const_cast<objectRegistry&>(mesh.thisDb());

    HashTable<GeoField*> flds(obr.lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

        bfld[patchi] == value;
    }
}


// ************************************************************************* //
