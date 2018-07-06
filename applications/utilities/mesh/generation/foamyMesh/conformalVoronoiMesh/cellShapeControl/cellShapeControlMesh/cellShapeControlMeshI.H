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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::cellShapeControlMesh::Vertex_handle Foam::cellShapeControlMesh::insert
(
    const Foam::point& pt,
    const scalar& size,
    const triad& alignment,
    const Foam::indexedVertexEnum::vertexType type
)
{
    Vertex_handle v = CellSizeDelaunay::insert
    (
        Point(pt.x(), pt.y(), pt.z())
    );
    v->type() = type;
    v->index() = getNewVertexIndex();
    v->procIndex() = Pstream::myProcNo();
    v->targetCellSize() = size;
    v->alignment() = tensor(alignment.x(), alignment.y(), alignment.z());

    return v;
}


Foam::cellShapeControlMesh::Vertex_handle Foam::cellShapeControlMesh::insertFar
(
    const Foam::point& pt
)
{
    Vertex_handle v = CellSizeDelaunay::insert
    (
        Point(pt.x(), pt.y(), pt.z())
    );
    v->type() = Vb::vtFar;
//    v->type() = Vb::vtExternalFeaturePoint;
    v->index() = getNewVertexIndex();
    v->procIndex() = Pstream::myProcNo();

    return v;
}


// ************************************************************************* //
