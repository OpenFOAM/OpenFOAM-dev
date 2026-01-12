/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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
    Foam::patchCutPlot

Description
    Functions for generating weights for a cut-plot of a patch

SourceFiles
    patchCutPlot.C

\*---------------------------------------------------------------------------*/

#ifndef patchCutPlot_H
#define patchCutPlot_H

#include "cutPlot.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class setWriter;

namespace patchCutPlot
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Structure containing the weight for a given cell and cut index
using cutPlot::weight;

//- Calculate and return weights
List<weight> calcWeights
(
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise = true
);

//- Calculate and return weights
template<class PatchType>
List<weight> calcWeights
(
    const PatchType& p,
    const scalarField& localPointXs,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise = true
);

//- Compute and return evenly-spaced cut coordinates
tmp<scalarField> calcCutXs
(
    const faceList& faces,
    const UList<vector>& faceAreas,
    const UList<vector>& faceNormals,
    const pointField& points,
    const scalarField& pointXs,
    const bool interpolate,
    const label nCuts,
    const label nIter,
    const bool debug = false,
    const word& functionName = word::null,
    const polyMesh& functionMesh = NullObjectRef<polyMesh>(),
    const setWriter& functionFormatter = NullObjectRef<setWriter>()
);

//- Compute and return evenly-spaced cut coordinates
template<class PatchType>
tmp<scalarField> calcCutXs
(
    const PatchType& p,
    const scalarField& localPointXs,
    const bool interpolate,
    const label nCuts,
    const label nIter,
    const bool debug = false,
    const word& functionName = word::null,
    const polyMesh& functionMesh = NullObjectRef<polyMesh>(),
    const setWriter& functionFormatter = NullObjectRef<setWriter>()
);

//- Write the layers as components of a tensor field for debugging
void writeLayers
(
    const SubList<face>& faces,
    const List<patchCutPlot::weight>& weights,
    const word& functionName,
    const fvMesh& functionMesh
);

//- Write the layers as components of a tensor field for debugging
void writeLayers
(
    const UIndirectList<face>& faces,
    const List<patchCutPlot::weight>& weights,
    const word& functionName,
    const fvMesh& functionMesh
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchCutPlot
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchCutPlotTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
