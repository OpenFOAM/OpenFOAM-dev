/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "polyMeshFilter.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class SetType>
void Foam::polyMeshFilter::updateSets(const mapPolyMesh& map)
{
    HashTable<const SetType*> sets =
        map.mesh().objectRegistry::lookupClass<const SetType>();

    forAllIter(typename HashTable<const SetType*>, sets, iter)
    {
        SetType& set = const_cast<SetType&>(*iter());
        set.updateMesh(map);
        set.sync(map.mesh());
    }

    IOobjectList Objects
    (
        map.mesh().time(),
        map.mesh().facesInstance(),
        "polyMesh/sets"
    );

    IOobjectList fileSets(Objects.lookupClass(SetType::typeName));

    forAllConstIter(IOobjectList, fileSets, iter)
    {
        if (!sets.found(iter.key()))
        {
            // Not in memory. Load it.
            SetType set(*iter());
            set.updateMesh(map);

            set.write();
        }
    }
}


template<class SetType>
void Foam::polyMeshFilter::copySets
(
    const polyMesh& oldMesh,
    const polyMesh& newMesh
)
{
    HashTable<const SetType*> sets =
        oldMesh.objectRegistry::lookupClass<const SetType>();

    forAllConstIter(typename HashTable<const SetType*>, sets, iter)
    {
        const SetType& set = *iter();

        if (newMesh.objectRegistry::foundObject<SetType>(set.name()))
        {
            const SetType& origSet =
                newMesh.objectRegistry::lookupObject<SetType>(set.name());

            const_cast<SetType&>(origSet) = set;
            const_cast<SetType&>(origSet).sync(newMesh);
        }
        else
        {
            SetType* newSet
            (
                new SetType(newMesh, set.name(), set, set.writeOpt())
            );

            newSet->store();
            newSet->sync(newMesh);
        }
    }
}


// ************************************************************************* //
