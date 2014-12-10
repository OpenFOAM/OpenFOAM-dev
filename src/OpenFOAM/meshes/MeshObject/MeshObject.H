/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::MeshObject

Description
    Templated abstract base-class for optional mesh objects used to automate
    their allocation to the mesh database and the mesh-modifier event-loop.

    MeshObject is templated on the type of mesh it is allocated to, the type of
    the mesh object (TopologicalMeshObject, GeometricMeshObject,
    MoveableMeshObject, UpdateableMeshObject) and the type of the actual object
    it is created for example:

    class leastSquaresVectors
    :
        public MeshObject<fvMesh, MoveableMeshObject, leastSquaresVectors>
    {
    .
    .
    .
        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
    };

    MeshObject types:

    TopologicalMeshObject: mesh object to be deleted on topology change
    GeometricMeshObject: mesh object to be deleted on geometry change
    MoveableMeshObject: mesh object to be updated in movePoints
    UpdateableMeshObject: mesh object to be updated in updateMesh or movePoints

    Note that movePoints must be provided for MeshObjects of type
    MoveableMeshObject and both movePoints and updateMesh functions must exist
    provided for MeshObjects of type UpdateableMeshObject.

SourceFiles
    MeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef MeshObject_H
#define MeshObject_H

#include "regIOobject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Class MeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh, template<class> class MeshObjectType, class Type>
class MeshObject
:
    public MeshObjectType<Mesh>
{

protected:

        // Reference to Mesh
        const Mesh& mesh_;


public:

    // Constructors

        explicit MeshObject(const Mesh& mesh);

        static const Type& New(const Mesh& mesh);

        template<class Data1>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1& d
        );

        template<class Data1, class Data2>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1&,
            const Data2&
        );

        template<class Data1, class Data2, class Data3>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1&,
            const Data2&,
            const Data3&
        );

        template<class Data1, class Data2, class Data3, class Data4>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1&,
            const Data2&,
            const Data3&,
            const Data4&
        );


    // Destructors

        virtual ~MeshObject();

        static bool Delete(const Mesh& mesh);


    // Member Functions

        const Mesh& mesh() const
        {
            return mesh_;
        }

        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        }
};


/*---------------------------------------------------------------------------*\
                           Class meshObject Declaration
\*---------------------------------------------------------------------------*/

class meshObject
:
    public regIOobject
{
public:

    // Declare name of the class and its debug switch
    ClassName("meshObject");

    // Constructors

        meshObject(const word& typeName, const objectRegistry& obr);


    // Static member functions

        template<class Mesh>
        static void movePoints(objectRegistry&);

        template<class Mesh>
        static void updateMesh(objectRegistry&, const mapPolyMesh&);

        template<class Mesh, template<class> class MeshObjectType>
        static void clear(objectRegistry&);

        //- Clear all meshObject derived from FromType up to (but not including)
        //  ToType. Used to clear e.g. all non-updateable meshObjects
        template
        <
            class Mesh,
            template<class> class FromType,
            template<class> class ToType
        >
        static void clearUpto(objectRegistry&);
};


/*---------------------------------------------------------------------------*\
                    Class TopologicalMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class TopologicalMeshObject
:
    public meshObject
{
public:

    TopologicalMeshObject(const word& typeName, const objectRegistry& obr)
    :
        meshObject(typeName, obr)
    {}
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

    GeometricMeshObject(const word& typeName, const objectRegistry& obr)
    :
        TopologicalMeshObject<Mesh>(typeName, obr)
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

    MoveableMeshObject(const word& typeName, const objectRegistry& obr)
    :
        GeometricMeshObject<Mesh>(typeName, obr)
    {}

    virtual bool movePoints() = 0;
};


/*---------------------------------------------------------------------------*\
                    Class UpdateableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class UpdateableMeshObject
:
    public MoveableMeshObject<Mesh>
{
public:

    UpdateableMeshObject(const word& typeName, const objectRegistry& obr)
    :
        MoveableMeshObject<Mesh>(typeName, obr)
    {}

    virtual void updateMesh(const mapPolyMesh& mpm) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeshObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
