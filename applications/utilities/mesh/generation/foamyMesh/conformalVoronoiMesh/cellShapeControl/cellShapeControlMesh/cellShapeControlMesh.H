/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::cellShapeControlMesh

Description

SourceFiles
    cellShapeControlMeshI.H
    cellShapeControlMesh.C

\*---------------------------------------------------------------------------*/

#ifndef cellShapeControlMesh_H
#define cellShapeControlMesh_H

#include "Time.H"
#include "scalar.H"
#include "point.H"
#include "tensor.H"
#include "triad.H"
#include "fileName.H"
#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "DistributedDelaunayMesh.H"
#include "CGALTriangulation3Ddefs.H"
#include "backgroundMeshDecomposition.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellSizeAndAlignmentControls;

/*---------------------------------------------------------------------------*\
                      Class cellShapeControlMesh Declaration
\*---------------------------------------------------------------------------*/

class cellShapeControlMesh
:
    public DistributedDelaunayMesh<CellSizeDelaunay>
{
public:

        typedef CellSizeDelaunay::Cell_handle      Cell_handle;
        typedef CellSizeDelaunay::Vertex_handle    Vertex_handle;
        typedef CellSizeDelaunay::Point            Point;


private:

        // Private data

        const Time& runTime_;

        mutable Cell_handle oldCellHandle_;

        const scalar defaultCellSize_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellShapeControlMesh(const cellShapeControlMesh&);

        //- Disallow default bitwise assignment
        void operator=(const cellShapeControlMesh&);


public:

    //- Runtime type information
    ClassName("cellShapeControlMesh");

    //- Return the mesh sub-directory name (usually "cellShapeControlMesh")
    static word meshSubDir;


    // Constructors

        explicit cellShapeControlMesh(const Time& runTime);


    //- Destructor
    ~cellShapeControlMesh();


    // Member Functions

        // Query

            const Time& time() const
            {
                return runTime_;
            }

            //- Calculate and return the barycentric coordinates for
            //  interpolating quantities on the background mesh
            void barycentricCoords
            (
                const Foam::point& pt,
                barycentric& bary,
                Cell_handle& ch
            ) const;

            boundBox bounds() const;


        // Edit

            label removePoints();

            //- Get the centres of all the tets
            tmp<pointField> cellCentres() const;

            inline Vertex_handle insert
            (
                const Foam::point& pt,
                const scalar& size,
                const triad& alignment,
                const Foam::indexedVertexEnum::vertexType type = Vb::vtInternal
            );

            inline Vertex_handle insertFar
            (
                const Foam::point& pt
            );

            void distribute(const backgroundMeshDecomposition& decomposition);

            tensorField dumpAlignments() const;

            void writeTriangulation();

            void write() const;

            label estimateCellCount
            (
                const autoPtr<backgroundMeshDecomposition>& decomposition
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellShapeControlMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
