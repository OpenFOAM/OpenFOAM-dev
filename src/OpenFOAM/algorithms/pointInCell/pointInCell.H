/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Function for determining if a point is within a cell of a polyMesh

SourceFiles
    pointInCell.C

\*---------------------------------------------------------------------------*/

#ifndef pointInCell_H
#define pointInCell_H

#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Enumeration for the sub-shapes used to perform the inside calculation
enum class pointInCellShapes
{
    facePlanes,
    faceCentreTris,
    faceDiagonalTris,
    tets
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Test if a point is in a given cell. For each of the cell's faces, define a
//  plane based on the face centre and the outward pointing face normal. The
//  point is considered to be inside the cell if it is below all such planes.
//  This is not guaranteed to be correct if the cell is not convex. Even if the
//  cell is convex, this approach is still somewhat approximate if the faces
//  are not flat.
bool pointInCellFacePlanes
(
    const point& p,
    const polyMesh& mesh,
    const label celli
);

//- Test if a point is in a given cell. For each of the cell's faces,
//  triangulate by connecting each edge to the face centre. Then proceed as for
//  pointInCellFacePlanes, but using these triangles instead of the faces.
//  Triangles are flat by definition, so this approach resolves the ambiguities
//  associated with representing faces as planes. It doesn't fix the convexity
//  requirement, however. If anything it makes it worse, as now the cell has
//  to be convex with respect to the triangles, of which there are more than
//  there are faces.
bool pointInCellFaceCentreTris
(
    const point& p,
    const polyMesh& mesh,
    const label celli
);

//- Test if a point is in a given cell. For each of the cell's faces,
//  triangulate by connecting each edge not containing the tetBasePt to the
//  tetBasePt. Then proceed as for pointInCellFacePlanes, but using these
//  triangles instead of the faces. This constructs fewer faces than
//  pointInCellFaceCentreTris, so it should be computationally advantageous
//  over that method. It shares the same limitations otherwise.
bool pointInCellFaceDiagTris
(
    const point& p,
    const polyMesh& mesh,
    const label celli
);

//- Test if a point is in a given cell by decomposing the cell into tetrahedra
//  and asking the tetrahedral decomposition engine to determine whether the
//  point is inside one of the tetrahedra. This method has the greatest cost,
//  but is the most general and should be robust to concave non-flat faces and
//  concave cells.
bool pointInCellTets
(
    const point& p,
    const polyMesh& mesh,
    const label celli
);

//- Test if a point is in a given cell
bool pointInCell
(
    const point& p,
    const polyMesh& mesh,
    const label celli,
    const pointInCellShapes = pointInCellShapes::tets
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
