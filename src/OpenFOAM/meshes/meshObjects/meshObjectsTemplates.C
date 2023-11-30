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

\*---------------------------------------------------------------------------*/

#include "meshObjects.H"
#include "MeshObjects.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Mesh>
void Foam::meshObjects::Delete(regIOobject& io)
{
    if (meshObjects::debug)
    {
        Pout<< "    Destroying " << io.name() << endl;
    }

    if (io.ownedByRegistry())
    {
        io.checkOut();
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to checkout and delete object " << io.name()
            << " not owned by registry."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Mesh>
void Foam::meshObjects::movePoints(objectRegistry& obr)
{
    HashTable<DeletableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<DeletableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::movePoints(objectRegistry&) :"
            << " moving " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<DeletableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<MoveableMeshObject<Mesh>>(*iter()))
        {
            if (meshObjects::debug)
            {
                Pout<< "    Moving " << iter()->io_.name() << endl;
            }
            dynamic_cast<MoveableMeshObject<Mesh>*>(iter())->movePoints();
        }
        else
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


template<class Mesh>
void Foam::meshObjects::distribute
(
    objectRegistry& obr,
    const polyDistributionMap& map
)
{
    HashTable<DeletableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<DeletableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::distribute(objectRegistry&, "
               "const polyDistributionMap& map) : updating "
            << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<DeletableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<TopoChangeableMeshObject<Mesh>>(*iter()))
        {
            if (meshObjects::debug)
            {
                Pout<< "    Distributing " << iter()->io_.name() << endl;
            }
            dynamic_cast<DistributeableMeshObject<Mesh>*>(iter())
            ->distribute(map);
        }
        else
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


template<class Mesh>
void Foam::meshObjects::topoChange
(
    objectRegistry& obr,
    const polyTopoChangeMap& map
)
{
    HashTable<DeletableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<DeletableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::topoChange(objectRegistry&, "
               "const polyTopoChangeMap& map) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<DeletableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<TopoChangeableMeshObject<Mesh>>(*iter()))
        {
            if (meshObjects::debug)
            {
                Pout<< "    Updating " << iter()->io_.name() << endl;
            }
            dynamic_cast<TopoChangeableMeshObject<Mesh>*>(iter())
                ->topoChange(map);
        }
        else
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


template<class Mesh>
void Foam::meshObjects::mapMesh
(
    objectRegistry& obr,
    const polyMeshMap& map
)
{
    HashTable<DeletableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<DeletableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::mapMesh(objectRegistry&, "
               "const polyMeshMap& map) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<DeletableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<TopoChangeableMeshObject<Mesh>>(*iter()))
        {
            if (meshObjects::debug)
            {
                Pout<< "    Updating " << iter()->io_.name() << endl;
            }
            dynamic_cast<TopoChangeableMeshObject<Mesh>*>(iter())->mapMesh(map);
        }
        else
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


template<class Mesh>
void Foam::meshObjects::addPatch(objectRegistry& obr, const label patchi)
{
    HashTable<DeletableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<DeletableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::addPatch(objectRegistry&, "
               "const label patchi) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<DeletableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<RepatchableMeshObject<Mesh>>(*iter()))
        {
            if (meshObjects::debug)
            {
                Pout<< "    Adding patch to " << iter()->io_.name() << endl;
            }
            dynamic_cast<RepatchableMeshObject<Mesh>*>(iter())
                ->addPatch(patchi);
        }
        else
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


template<class Mesh>
void Foam::meshObjects::reorderPatches
(
    objectRegistry& obr,
    const labelUList& newToOld,
    const bool validBoundary
)
{
    HashTable<DeletableMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<DeletableMeshObject<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::addPatch(objectRegistry&, "
               "const labelUList&, const bool) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<DeletableMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<RepatchableMeshObject<Mesh>>(*iter()))
        {
            if (meshObjects::debug)
            {
                Pout<< "    Adding patch to " << iter()->io_.name() << endl;
            }
            dynamic_cast<RepatchableMeshObject<Mesh>*>(iter())->reorderPatches
            (
                newToOld,
                validBoundary
            );
        }
        else
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


template<class Mesh, template<class> class MeshObjectType>
void Foam::meshObjects::clear(objectRegistry& obr)
{
    HashTable<MeshObjectType<Mesh>*> meshObjects
    (
        obr.lookupClass<MeshObjectType<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::clear(objectRegistry&) :"
            << " clearing " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter(typename HashTable<MeshObjectType<Mesh>*>, meshObjects, iter)
    {
        Delete<Mesh>(iter()->io_);
    }
}


template
<
    class Mesh,
    template<class> class FromType,
    template<class> class ToType
>
void Foam::meshObjects::clearUpto(objectRegistry& obr)
{
    HashTable<FromType<Mesh>*> meshObjects
    (
        obr.lookupClass<FromType<Mesh>>()
    );

    if (meshObjects::debug)
    {
        Pout<< "meshObjects::clearUpto(objectRegistry&) :"
            << " clearing " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter(typename HashTable<FromType<Mesh>*>, meshObjects, iter)
    {
        if (!isA<ToType<Mesh>>(*iter()))
        {
            Delete<Mesh>(iter()->io_);
        }
    }
}


// ************************************************************************* //
