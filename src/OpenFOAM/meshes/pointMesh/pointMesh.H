/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "MeshObject.H"
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
    public MeshObject<polyMesh, PatchMeshObject, pointMesh>,
    public GeoMesh<polyMesh>
{
    // Permanent data

        //- Boundary mesh
        pointBoundaryMesh boundary_;


    // Private Member Functions

        //- Map all fields
        void mapFields(const mapPolyMesh&);


public:

    // Declare name of the class and its debug switch
    ClassName("pointMesh");


    typedef pointMesh Mesh;
    typedef pointBoundaryMesh BoundaryMesh;


    // Constructors

        //- Construct from polyMesh
        explicit pointMesh(const polyMesh& pMesh);

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

        //- Reset for changed polyMesh
        void reset(const bool validBoundary);


        // Mesh callbacks

            //- Move points
            virtual bool movePoints();

            //- Update the mesh corresponding to given map
            virtual void updateMesh(const mapPolyMesh& mpm);

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
