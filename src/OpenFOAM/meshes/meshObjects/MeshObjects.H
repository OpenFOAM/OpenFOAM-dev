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

Class
    Foam::TopologicalMeshObject

Description
    MeshObject types:

    - TopologicalMeshObject: mesh object to be deleted on topology change
    - GeometricMeshObject: mesh object to be deleted on geometry change
    - MoveableMeshObject: mesh object to be updated in movePoints
    - UpdateableMeshObject: mesh object to be updated in topoChange or
        movePoints
    - PatchMeshObject: mesh object to be additionally updated patch changes

\*---------------------------------------------------------------------------*/

#ifndef MeshObjects_H
#define MeshObjects_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;
class meshObjects;

/*---------------------------------------------------------------------------*\
                    Class TopologicalMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class TopologicalMeshObject
{
    //- Reference to the regIOobject of the base-class
    regIOobject& io_;

    //- The regIOobject reference is used by the meshObjects functions
    friend class meshObjects;

public:

    TopologicalMeshObject(regIOobject& io, const Mesh& mesh)
    :
        io_(io)
    {}

    //- Virtual destructor to make class polymorphic
    virtual ~TopologicalMeshObject() = default;
};


/*---------------------------------------------------------------------------*\
                     Class GeometricMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class GeometricMeshObject
:
    public TopologicalMeshObject<Mesh>
{
public:

    GeometricMeshObject(regIOobject& io, const Mesh& mesh)
    :
        TopologicalMeshObject<Mesh>(io, mesh)
    {}
};


/*---------------------------------------------------------------------------*\
                     Class MoveableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class MoveableMeshObject
:
    public GeometricMeshObject<Mesh>
{
public:

    MoveableMeshObject(regIOobject& io, const Mesh& mesh)
    :
        GeometricMeshObject<Mesh>(io, mesh)
    {}

    //- Update for mesh motion
    virtual bool movePoints() = 0;
};


/*---------------------------------------------------------------------------*\
                   Class DistributeableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class DistributeableMeshObject
:
    public MoveableMeshObject<Mesh>
{
public:

    DistributeableMeshObject(regIOobject& io, const Mesh& mesh)
    :
        MoveableMeshObject<Mesh>(io, mesh)
    {}

    //- Redistribute or update using the given distribution map
    virtual void distribute(const polyDistributionMap& map) = 0;
};


/*---------------------------------------------------------------------------*\
                    Class UpdateableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class UpdateableMeshObject
:
    public DistributeableMeshObject<Mesh>
{
public:

    UpdateableMeshObject(regIOobject& io, const Mesh& mesh)
    :
        DistributeableMeshObject<Mesh>(io, mesh)
    {}

    //- Update topology using the given map
    virtual void topoChange(const polyTopoChangeMap& map) = 0;

    //- Update from another mesh using the given map
    virtual void mapMesh(const polyMeshMap& map) = 0;
};


/*---------------------------------------------------------------------------*\
                       Class PatchMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class PatchMeshObject
:
    public UpdateableMeshObject<Mesh>
{
public:

    PatchMeshObject(regIOobject& io, const Mesh& mesh)
    :
        UpdateableMeshObject<Mesh>(io, mesh)
    {}

    //- Reordered/removed trailing patches. If validBoundary call is parallel
    //  synced and all add the same patch with same settings
    virtual void reorderPatches
    (
        const labelUList& newToOld,
        const bool validBoundary
    ) = 0;

    //- Inserted patch at patchi
    virtual void addPatch(const label patchi) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
