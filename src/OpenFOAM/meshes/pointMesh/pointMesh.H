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

Class
    Foam::pointMesh

Description
    Mesh representing a set of points created from polyMesh.

\*---------------------------------------------------------------------------*/

#ifndef pointMesh_H
#define pointMesh_H

#include "GeoMesh.H"
#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pointMesh Declaration
\*---------------------------------------------------------------------------*/

class pointMesh
:
    public DemandDrivenMeshObject<polyMesh, PatchMeshObject, pointMesh>,
    public GeoMesh<polyMesh>
{
    // Permanent data

        //- Boundary mesh
        pointBoundaryMesh boundary_;


protected:

    friend class DemandDrivenMeshObject<polyMesh, PatchMeshObject, pointMesh>;

    // Protected Constructors

        //- Construct from polyMesh
        explicit pointMesh(const polyMesh& pMesh);


public:

    // Public Typedefs

        typedef pointMesh Mesh;
        typedef pointBoundaryMesh BoundaryMesh;


    // Declare name of the class and its debug switch
    ClassName("pointMesh");


    // Static data

        //- Set of names of registered geometric fields
        const static HashSet<word> geometryFields;


    // Constructors

        //- Disallow default bitwise copy construction
        pointMesh(const pointMesh&) = delete;


    //- Destructor
    ~pointMesh();


    // Member Functions

        //- Return number of points
        label size() const
        {
            return size(*this);
        }

        //- Return number of points
        static label size(const Mesh& mesh)
        {
            return mesh.GeoMesh<polyMesh>::mesh_.nPoints();
        }

        //- Return reference to boundary mesh
        const pointBoundaryMesh& boundary() const
        {
            return boundary_;
        }

        //- Return parallel info
        const globalMeshData& globalData() const
        {
            return GeoMesh<polyMesh>::mesh_.globalData();
        }

        //- Return database. For now is its polyMesh.
        const objectRegistry& thisDb() const
        {
            return GeoMesh<polyMesh>::mesh_.thisDb();
        }

        //- Reset pointMesh with respect to the updated polyMesh
        //  For run-time mesh replacement and mesh to mesh mapping
        void reset();


        // Mesh callbacks

            //- Move points
            virtual bool movePoints();

            //- Update the mesh corresponding to given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Reordered/removed trailing patches. If validBoundary call is
            //  parallel synced and all add the same patch with same settings
            virtual void reorderPatches
            (
                const labelUList& newToOld,
                const bool validBoundary
            );

            //- Inserted patch at patchi
            virtual void addPatch(const label patchi);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointMesh&) = delete;

        bool operator!=(const pointMesh& pm) const
        {
            return &pm != this;
        }

        bool operator==(const pointMesh& pm) const
        {
            return &pm == this;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
