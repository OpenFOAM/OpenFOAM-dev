/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::DelaunayMeshTools

Description
    Collection of functions for operating on a Delaunay mesh. Includes:

      - Functions for writing to an OBJ file
      - Functions for extracting fields from the Delaunay triangulation

SourceFiles
    DelaunayMeshToolsI.H
    DelaunayMeshTools.C

\*---------------------------------------------------------------------------*/

#ifndef DelaunayMeshTools_H
#define DelaunayMeshTools_H

#include "fileName.H"
#include "List.H"
#include "point.H"
#include "CGALTriangulation3Ddefs.H"
#include "indexedVertexEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Namespace DelaunayMeshTools Declaration
\*---------------------------------------------------------------------------*/

namespace DelaunayMeshTools
{

// OBJ writing

    //- Write list of points to file
    void writeOBJ(const fileName& fName, const List<Foam::point>& points);

    //- Write list of points to file
    void writeOBJ(const fileName& fName, const List<Vb>& points);

    //- Write an OBJ mesh consisting of points and faces
    void writeObjMesh
    (
        const fileName& fName,
        const pointField& points,
        const faceList& faces
    );

    //- Write Delaunay points in the range between (and including)
    //  type startPointType and endPointType to an OBJ file
    template<class Triangulation>
    void writeOBJ
    (
        const fileName& fName,
        const Triangulation& t,
        const indexedVertexEnum::vertexType startPointType,
        const indexedVertexEnum::vertexType endPointType
    );

    //- Write Delaunay points of type pointType to .obj file
    template<class Triangulation>
    void writeOBJ
    (
        const fileName& fName,
        const Triangulation& t,
        const indexedVertexEnum::vertexType pointType
    );

    //- Write the fixed Delaunay points to an OBJ file
    template<class Triangulation>
    void writeFixedPoints(const fileName& fName, const Triangulation& t);

    //- Write the boundary Delaunay points to an OBJ file
    template<class Triangulation>
    void writeBoundaryPoints(const fileName& fName, const Triangulation& t);

    //- Write the processor interface to an OBJ file
    template<class Triangulation>
    void writeProcessorInterface
    (
        const fileName& fName,
        const Triangulation& t,
        const faceList& faces
    );

    //- Write the internal Delaunay vertices of the tessellation as a
    //  pointField that may be used to restart the meshing process
    template<class Triangulation>
    void writeInternalDelaunayVertices
    (
        const fileName& instance,
        const Triangulation& t
    );

    //- Draws a tet cell to an output stream. The offset is supplied as the tet
    //  number to be drawn.
    template<class CellHandle>
    void drawDelaunayCell(Ostream& os, const CellHandle& c, label offset = 0);


// Field extraction

    //- Extract all points in vertex-index order
    template<class Triangulation>
    tmp<pointField> allPoints(const Triangulation& t);


} // End namespace DelaunayMeshTools


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DelaunayMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
