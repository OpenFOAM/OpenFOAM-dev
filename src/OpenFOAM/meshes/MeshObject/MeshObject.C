/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "MeshObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::MeshObject(const Mesh& mesh)
:
    MeshObjectType<Mesh>(Type::typeName, mesh.thisDb()),
    mesh_(mesh)
{}


template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::MeshObject
(
    const Mesh& mesh,
    const IOobject& io
)
:
    MeshObjectType<Mesh>(Type::typeName, io),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    Mesh& mesh
)
{
    if
    (
        mesh.thisDb().objectRegistry::template foundObject<Type>
        (
            Type::typeName
        )
    )
    {
        return mesh.thisDb().objectRegistry::template lookupObjectRef<Type>
        (
            Type::typeName
        );
    }
    else
    {
        if (meshObject::debug)
        {
            Pout<< "MeshObject::New(" << Mesh::typeName
                << "&) : constructing " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(mesh);

        regIOobject::store(static_cast<MeshObjectType<Mesh>*>(objectPtr));

        return *objectPtr;
    }
}


template<class Mesh, template<class> class MeshObjectType, class Type>
const Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    const Mesh& mesh
)
{
    if
    (
        mesh.thisDb().objectRegistry::template foundObject<Type>
        (
            Type::typeName
        )
    )
    {
        return mesh.thisDb().objectRegistry::template lookupObjectRef<Type>
        (
            Type::typeName
        );
    }
    else
    {
        if (meshObject::debug)
        {
            Pout<< "MeshObject::New(const " << Mesh::typeName
                << "&) : constructing " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(mesh);

        regIOobject::store(static_cast<MeshObjectType<Mesh>*>(objectPtr));

        return *objectPtr;
    }
}


template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    Mesh& mesh,
    const Args&... args
)
{
    if
    (
        mesh.thisDb().objectRegistry::template foundObject<Type>
        (
            Type::typeName
        )
    )
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
    else
    {
        if (meshObject::debug)
        {
           Pout<< "MeshObject::New(" << Mesh::typeName
                << "&, const Data1&) : constructing " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(mesh, args...);

        regIOobject::store(static_cast<MeshObjectType<Mesh>*>(objectPtr));

        return *objectPtr;
    }
}


template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
const Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    const Mesh& mesh,
    const Args&... args
)
{
    if
    (
        mesh.thisDb().objectRegistry::template foundObject<Type>
        (
            Type::typeName
        )
    )
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
    else
    {
        if (meshObject::debug)
        {
           Pout<< "MeshObject::New(const " << Mesh::typeName
                << "&, const Data1&) : constructing " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(mesh, args...);

        regIOobject::store(static_cast<MeshObjectType<Mesh>*>(objectPtr));

        return *objectPtr;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
bool Foam::MeshObject<Mesh, MeshObjectType, Type>::Delete(const Mesh& mesh)
{
    if
    (
        mesh.thisDb().objectRegistry::template foundObject<Type>
        (
            Type::typeName
        )
    )
    {
        if (meshObject::debug)
        {
            Pout<< "MeshObject::Delete(const Mesh&) : deleting "
                << Type::typeName << endl;
        }

        return mesh.thisDb().checkOut
        (
            const_cast<Type&>
            (
                mesh.thisDb().objectRegistry::template lookupObject<Type>
                (
                    Type::typeName
                )
            )
        );
    }
    else
    {
        return false;
    }
}


template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::~MeshObject()
{
    MeshObjectType<Mesh>::release();
}


template<class Mesh>
void Foam::meshObject::movePoints(objectRegistry& obr)
{
    HashTable<GeometricMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::movePoints(objectRegistry&) :"
            << " moving " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<GeometricMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<MoveableMeshObject<Mesh>>(*iter()))
        {
            if (meshObject::debug)
            {
                Pout<< "    Moving " << iter()->name() << endl;
            }
            dynamic_cast<MoveableMeshObject<Mesh>*>(iter())->movePoints();
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << iter()->name() << endl;
            }
            obr.checkOut(*iter());
        }
    }
}


template<class Mesh>
void Foam::meshObject::updateMesh(objectRegistry& obr, const mapPolyMesh& mpm)
{
    HashTable<GeometricMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::updateMesh(objectRegistry&, "
               "const mapPolyMesh& mpm) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<GeometricMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<UpdateableMeshObject<Mesh>>(*iter()))
        {
            if (meshObject::debug)
            {
                Pout<< "    Updating " << iter()->name() << endl;
            }
            dynamic_cast<UpdateableMeshObject<Mesh>*>(iter())->updateMesh(mpm);
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << iter()->name() << endl;
            }
            obr.checkOut(*iter());
        }
    }
}


template<class Mesh>
void Foam::meshObject::addPatch(objectRegistry& obr, const label patchi)
{
    HashTable<GeometricMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::addPatch(objectRegistry&, "
               "const label patchi) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<GeometricMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<PatchMeshObject<Mesh>>(*iter()))
        {
            if (meshObject::debug)
            {
                Pout<< "    Adding patch to " << iter()->name() << endl;
            }
            dynamic_cast<PatchMeshObject<Mesh>*>(iter())->addPatch(patchi);
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << iter()->name() << endl;
            }
            obr.checkOut(*iter());
        }
    }
}


template<class Mesh>
void Foam::meshObject::reorderPatches
(
    objectRegistry& obr,
    const labelUList& newToOld,
    const bool validBoundary
)
{
    HashTable<GeometricMeshObject<Mesh>*> meshObjects
    (
        obr.lookupClass<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::addPatch(objectRegistry&, "
               "const labelUList&, const bool) : updating " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter
    (
        typename HashTable<GeometricMeshObject<Mesh>*>,
        meshObjects,
        iter
    )
    {
        if (isA<PatchMeshObject<Mesh>>(*iter()))
        {
            if (meshObject::debug)
            {
                Pout<< "    Adding patch to " << iter()->name() << endl;
            }
            dynamic_cast<PatchMeshObject<Mesh>*>(iter())->reorderPatches
            (
                newToOld,
                validBoundary
            );
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << iter()->name() << endl;
            }
            obr.checkOut(*iter());
        }
    }
}


template<class Mesh, template<class> class MeshObjectType>
void Foam::meshObject::clear(objectRegistry& obr)
{
    HashTable<MeshObjectType<Mesh>*> meshObjects
    (
        obr.lookupClass<MeshObjectType<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::clear(objectRegistry&) :"
            << " clearing " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter(typename HashTable<MeshObjectType<Mesh>*>, meshObjects, iter)
    {
        if (meshObject::debug)
        {
            Pout<< "    Destroying " << iter()->name() << endl;
        }
        obr.checkOut(*iter());
    }
}


template
<
    class Mesh,
    template<class> class FromType,
    template<class> class ToType
>
void Foam::meshObject::clearUpto(objectRegistry& obr)
{
    HashTable<FromType<Mesh>*> meshObjects
    (
        obr.lookupClass<FromType<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::clearUpto(objectRegistry&) :"
            << " clearing " << Mesh::typeName
            << " meshObjects for region " << obr.name() << endl;
    }

    forAllIter(typename HashTable<FromType<Mesh>*>, meshObjects, iter)
    {
        if (!isA<ToType<Mesh>>(*iter()))
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << iter()->name() << endl;
            }
            obr.checkOut(*iter());
        }
    }
}


// ************************************************************************* //
