/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "DemandDrivenMeshObject.H"
#include "meshObjects.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
DemandDrivenMeshObject
(
    const IOobject& io,
    const Mesh& mesh
)
:
    RegIOobject(io),
    MeshObjectType<Mesh>(*this),
    mesh_(mesh)
{}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
DemandDrivenMeshObject
(
    const word& name,
    const Mesh& mesh
)
:
    RegIOobject
    (
        IOobject
        (
            name,
            mesh.db().instance(),
            mesh.db()
        )
    ),
    MeshObjectType<Mesh>(*this),
    mesh_(mesh)
{}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
DemandDrivenMeshObject
(
    const Mesh& mesh
)
:
    DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>
    (
        Type::typeName,
        mesh
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
Type& Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::New
(
    const word& name,
    const Mesh& mesh
)
{
    if (found(name, mesh))
    {
        return mesh.db().objectRegistry::template lookupObjectRef<Type>
        (
            name
        );
    }
    else
    {
        if (meshObjects::debug)
        {
            Pout<< "DemandDrivenMeshObject::New(" << Mesh::typeName
                << "&) : constructing " << name
                << " of type " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(name, mesh);

        return regIOobject::store(objectPtr);
    }
}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
Type& Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::New
(
    const Mesh& mesh
)
{
    if (found(mesh))
    {
        return mesh.db().objectRegistry::template lookupObjectRef<Type>
        (
            Type::typeName
        );
    }
    else
    {
        if (meshObjects::debug)
        {
            Pout<< "DemandDrivenMeshObject::New(" << Mesh::typeName
                << "&) : constructing " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(mesh);

        return regIOobject::store(objectPtr);
    }
}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
template<class... Args>
Type& Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::New
(
    const word& name,
    const Mesh& mesh,
    const Args&... args
)
{
    if (found(name, mesh))
    {
        return mesh.db().objectRegistry::template lookupObjectRef<Type>
        (
            name
        );
    }
    else
    {
        if (meshObjects::debug)
        {
            Pout<< "DemandDrivenMeshObject::New(" << Mesh::typeName
                << "&, const Data1&) : constructing " << name
                << " of type " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(name, mesh, args...);

        return regIOobject::store(objectPtr);
    }
}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
template<class... Args>
Type& Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::New
(
    const Mesh& mesh,
    const Args&... args
)
{
    if (found(mesh))
    {
        return mesh.db().objectRegistry::template lookupObjectRef<Type>
        (
            Type::typeName
        );
    }
    else
    {
        if (meshObjects::debug)
        {
            Pout<< "DemandDrivenMeshObject::New(" << Mesh::typeName
                << "&, const Data1&) : constructing " << Type::typeName
                << " for region " << mesh.name() << endl;
        }

        Type* objectPtr = new Type(mesh, args...);

        return regIOobject::store(objectPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
~DemandDrivenMeshObject()
{
    regIOobject::release();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
const Foam::word&
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
type() const
{
    return Type::typeName;
}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
bool Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
found
(
    const word& name,
    const Mesh& mesh
)
{
    return mesh.db().objectRegistry::template foundObject<Type>(name);
}


template
<
    class Mesh,
    template<class> class MeshObjectType,
    class Type,
    class RegIOobject
>
bool Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type, RegIOobject>::
found
(
    const Mesh& mesh
)
{
    return found(Type::typeName, mesh);
}


// ************************************************************************* //
