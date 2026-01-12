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
    Foam::cellCutPlot

Description
    Functions for generating weights for a cut-plot of a patch

SourceFiles
    cellCutPlot.C

\*---------------------------------------------------------------------------*/

#ifndef cellCutPlot_H
#define cellCutPlot_H

#include "cutPlot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class generatedCellZone;
class setWriter;

namespace cellCutPlot
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Structure containing the weight for a given cell and cut index
using cutPlot::weight;

//- Calculate and return weights
List<weight> calcWeights
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise = true
);

//- Compute and return evenly-spaced cut coordinates
tmp<scalarField> calcCutXs
(
    const polyMesh& mesh,
    const generatedCellZone& zone,
    const scalarField& pointXs,
    const bool interpolate,
    const label nCuts,
    const label nIter,
    const bool debug = false,
    const word& functionName = word::null,
    const setWriter& functionFormatter = NullObjectRef<setWriter>()
);

//- Write the layers as components of a tensor field for debugging
void writeLayers
(
    const fvMesh& mesh,
    const List<cellCutPlot::weight>& weights,
    const word& functionName
);

//- Construct plot values from cell values given a set of weights
template<class Type>
tmp<Field<Type>> applyWeights
(
    const label n,
    const List<cellCutPlot::weight>& weights,
    const Field<Type>& cellValues
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCutPlot
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
