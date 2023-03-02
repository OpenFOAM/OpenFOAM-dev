/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

Namespace
    Foam::cutPoly

Description
    Low level functions for cutting poly faces and cells

    Cutting a face using an iso-surface defined at the points is done by
    walking around the face. A start point is chosen and then the walk proceeds
    until a "cut" edge is found with points on either side of the iso-surface.
    This edge is noted and the walk continues until another cut edge is found.
    These two edge cuts, along with all points in-between, form a sub-face on
    one side of the iso-surface. This algorithm continues around the entire
    face, potentially cutting off multiple sub-faces. The remainder is a single
    contiguous face on the other side of the iso-surface.

    Cutting a cell is then done by walking from cut-edge to cut-edge around the
    cell's faces.

SourceFiles
    cutPoly.C

\*---------------------------------------------------------------------------*/

#ifndef cutPoly_H
#define cutPoly_H

#include "cellEdgeAddressing.H"
#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OBJstream;

namespace cutPoly
{

//- This flag determines which side of the iso-surface is separated into
//  multiple sub-faces, and which side is kept as the central contiguous face
static const bool separatedBelow = true;

//- Return the cuts for a given face. This returns a list of pairs of
//  face-edge indices. Both edges in the pair are cut, and connecting these
//  cuts creates an edge of the iso-surface.
List<labelPair> faceCuts
(
    const face& f,
    const scalarField& pAlphas,
    const scalar isoAlpha
);

//- Return the cuts for a given cell. This returns a list of lists of cell-edge
//  indices. Each sub-list represents a single closed loop of cell-edges. The
//  edges listed are all cut, and connecting these cuts up creates a face of
//  the iso-surface.
labelListList cellCuts
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const faceUList& fs,
    const List<List<labelPair>>& fCuts,
    const scalarField& pAlphas,
    const scalar isoAlpha
);

//- Write the cuts for a given face to an OBJ file
void writeFaceCuts
(
    const face& f,
    const List<labelPair>& fCuts,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    OBJstream& obj
);

//- Write the cuts for a given cell to an OBJ file
void writeCellCuts
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const List<List<label>>& cCuts,
    const faceList& fs,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    OBJstream& obj
);

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
