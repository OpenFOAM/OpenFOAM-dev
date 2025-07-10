/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::meshBoundarySearch

Description
    Mesh object that implements searches within the local boundary faces

SourceFiles
    meshBoundarySearch.C

\*---------------------------------------------------------------------------*/

#ifndef meshBoundarySearch_H
#define meshBoundarySearch_H

#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "treeDataFace.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshBoundarySearch Declaration
\*---------------------------------------------------------------------------*/

class meshBoundarySearch
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        DeletableMeshObject,
        meshBoundarySearch
    >
{
    // Private Classes

        // Intersection operation class that checks previous successful hits so
        // that they are not duplicated
        class findUniqueIntersectOp;


    // Private Data

        //- Boundary tree
        const indexedOctree<treeDataFace> boundaryTree_;


public:

    //- Run-time type information
    TypeName("meshBoundarySearch");


    // Constructors

        //- Construct from mesh
        meshBoundarySearch(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        meshBoundarySearch(const meshBoundarySearch&) = delete;


    //- Destructor
    virtual ~meshBoundarySearch();


    // Member Functions

        // Access

            //- Access the boundary tree
            inline const indexedOctree<treeDataFace>& boundaryTree() const
            {
                return boundaryTree_;
            }


        // Queries

            //- Find the boundary face with centre closest to the given point
            label findNearestBoundaryFace(const point& p) const;

            //- Find the first intersection of the boundary with the line
            //  defined by the given start and end points
            pointIndexHit intersection
            (
                const point& pStart,
                const point& pEnd
            ) const;

            //- Find all intersections of the boundary with the line
            //  defined by the given start and end points
            List<pointIndexHit> intersections
            (
                const point& pStart,
                const point& pEnd
            ) const;

            //- Determine whether the given point is inside the mesh boundary
            bool isInside(const point&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshBoundarySearch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
