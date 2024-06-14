/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::DemandDrivenMeshObject

Description
    Templated abstract base-class for demand-driven mesh objects used to
    automate their allocation to the mesh database and the mesh-modifier
    event-loop.

    DemandDrivenMeshObject is templated on the type of mesh it is allocated to,
    the type of the mesh object (DeletableMeshObject, MoveableMeshObject,
    DistributeableMeshObject, TopoChangeableMeshObject) and the type of the
    actual object it is created for example:

    \verbatim
    class leastSquaresVectors
    :
        public DemandDrivenMeshObject
        <
            fvMesh,
            MoveableMeshObject,
            leastSquaresVectors
        >
    {
    .
    .
    .
        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
    };
    \endverbatim

    MeshObject types:

    - DeletableMeshObject:
        mesh object to be deleted after any mesh change
    - MoveableMeshObject:
        mesh object to be updated after mesh motion otherwise deleted
    - DistributeableMeshObject
        mesh object to be updated after mesh redistribution or motion
        otherwise deleted
    - TopoChangeableMeshObject:
        mesh object to be updated after mesh topology change,
        mesh-to-mesh mapping, redistribution or motion otherwise deleted
    - RepatchableMeshObject:
        mesh object to be updated on patch or topology change,
        mesh-to-mesh mapping, redistribution or motion otherwise deleted

    DemandDrivenMeshObject should always be constructed and accessed via the New
    methods provided so that they are held and maintained by the objectRegistry.
    To ensure this use constructors of the concrete derived types should be
    private or protected and friendship with the DemandDrivenMeshObject
    base-class declared so that the New functions can call the the constructors.

SourceFiles
    DemandDrivenMeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef DemandDrivenMeshObject_H
#define DemandDrivenMeshObject_H

#include "MeshObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class DemandDrivenMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh, template<class> class MeshObjectType, class Type>
class DemandDrivenMeshObject
:
    public regIOobject,
    public MeshObjectType<Mesh>
{
    // Private member data

        // Reference to Mesh
        const Mesh& mesh_;


protected:

    // Constructors

        //- Construct from mesh and IOobject
        //  Only derived classes can construct DemandDrivenMeshObject
        DemandDrivenMeshObject(const IOobject& io, const Mesh& mesh);

        //- Construct from mesh and name
        //  Only derived classes can construct DemandDrivenMeshObject
        DemandDrivenMeshObject(const word& name, const Mesh& mesh);

        //- Construct from mesh, the name is set to Type::typeName
        //  Only derived classes can construct DemandDrivenMeshObject
        DemandDrivenMeshObject(const Mesh& mesh);


public:

    //- Runtime type information
    virtual const word& type() const;


    // Constructors

        //- Construct and return the named DemandDrivenMeshObject
        static Type& New(const word& name, const Mesh& mesh);

        //- Construct and return the DemandDrivenMeshObject named Type::typeName
        static Type& New(const Mesh& mesh);

        //- Construct and return the named DemandDrivenMeshObject
        //  with the additional arguments
        template<class... Args>
        static Type& New
        (
            const word& name,
            const Mesh& mesh,
            const Args&... args
        );

        //- Construct and return the DemandDrivenMeshObject named Type::typeName
        //  with the additional arguments
        template<class... Args>
        static Type& New
        (
            const Mesh& mesh,
            const Args&... args
        );


    // Destructors

        virtual ~DemandDrivenMeshObject();


    // Member Functions

        //- Return true if the DemandDrivenMeshObject with the given name
        //  is found in the mesh registry
        static bool found(const word& name, const Mesh& mesh);

        //- Return true if the DemandDrivenMeshObject named Type::typeName
        //  is found in the mesh registry
        static bool found(const Mesh& mesh);

        const Mesh& mesh() const
        {
            return mesh_;
        }

        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DemandDrivenMeshObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
