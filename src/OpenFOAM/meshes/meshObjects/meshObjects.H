/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::meshObjects

Description
    Templated abstract base-class for optional mesh objects used to automate
    their allocation to the mesh database and the mesh-modifier event-loop.

    MeshObject is templated on the type of mesh it is allocated to, the type of
    the mesh object (TopologicalMeshObject, GeometricMeshObject,
    MoveableMeshObject, DistributeableMeshObject, UpdateableMeshObject) and the
    type of the actual object it is created for example:

    \verbatim
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
    \endverbatim

    MeshObject types:

    - TopologicalMeshObject: mesh object to be deleted on topology change
    - GeometricMeshObject: mesh object to be deleted on geometry change
    - MoveableMeshObject: mesh object to be updated in movePoints
    - UpdateableMeshObject: mesh object to be updated in topoChange or
        movePoints
    - PatchMeshObject: mesh object to be additionally updated patch changes

    Note:
        movePoints must be provided for MeshObjects of type MoveableMeshObject
        and both movePoints and topoChange functions must exist, provided for
        MeshObjects of type UpdateableMeshObject.

SourceFiles
    meshObject.C
    meshObjects.C

\*---------------------------------------------------------------------------*/

#ifndef meshObjects_H
#define meshObjects_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                          Class meshObjects Declaration
\*---------------------------------------------------------------------------*/

class meshObjects
{
private:

    //- Checkout and delete object if owned by registry
    //  otherwise error
    template<class Mesh>
    static void Delete(regIOobject&);


public:

    ClassName("meshObjects");

    template<class Mesh>
    static void movePoints(objectRegistry&);

    template<class Mesh>
    static void topoChange(objectRegistry&, const polyTopoChangeMap&);

    template<class Mesh>
    static void mapMesh(objectRegistry&, const polyMeshMap&);

    template<class Mesh>
    static void distribute(objectRegistry&, const polyDistributionMap&);

    template<class Mesh>
    static void addPatch(objectRegistry&, const label patchi);

    template<class Mesh>
    static void reorderPatches
    (
        objectRegistry&,
        const labelUList& newToOld,
        const bool validBoundary
    );

    template<class Mesh, template<class> class MeshObjectType>
    static void clear(objectRegistry&);

    //- Clear all meshObjects derived from FromType up to (but not including)
    //  ToType. Used to clear e.g. all non-updateable meshObjects
    template
    <
        class Mesh,
        template<class> class FromType,
        template<class> class ToType
    >
    static void clearUpto(objectRegistry&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshObjectsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
