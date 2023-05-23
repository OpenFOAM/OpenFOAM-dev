/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::DemandDrivenMeshObject
(
    const Mesh& mesh
)
:
    regIOobject
    (
        IOobject
        (
            Type::typeName,
            mesh.thisDb().instance(),
            mesh.thisDb()
        )
    ),
    MeshObjectType<Mesh>(*this, mesh),
    mesh_(mesh)
{}


template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::DemandDrivenMeshObject
(
    const Mesh& mesh,
    const IOobject& io
)
:
    regIOobject(io),
    MeshObjectType<Mesh>(*this, mesh),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
Type& Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::New
(
    const Mesh& mesh
)
{
    if (found(mesh))
    {
        return mesh.thisDb().objectRegistry::template lookupObjectRef<Type>
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


template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
Type& Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::New
(
    const Mesh& mesh,
    const Args&... args
)
{
    if (found(mesh))
    {
        return mesh.thisDb().objectRegistry::template lookupObjectRef<Type>
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

template<class Mesh, template<class> class MeshObjectType, class Type>
bool Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::Delete
(
    const Mesh& mesh
)
{
    if (found(mesh))
    {
        if (meshObjects::debug)
        {
            Pout<< "DemandDrivenMeshObject::Delete(const Mesh&) : deleting "
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
Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::
~DemandDrivenMeshObject()
{
    regIOobject::release();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
bool Foam::DemandDrivenMeshObject<Mesh, MeshObjectType, Type>::found
(
    const Mesh& mesh
)
{
    return mesh.thisDb().objectRegistry::template foundObject<Type>
    (
        Type::typeName
    );
}


// ************************************************************************* //
