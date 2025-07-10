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

\*---------------------------------------------------------------------------*/

#include "lineCellFace.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(lineCellFace, 0);
    addToRunTimeSelectionTable(sampledSet, lineCellFace, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::lineCellFace::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>& samplingDistances,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    lineFace::calcSamples
    (
        mesh(),
        start_,
        end_,
        storeBothProcFaces_ ? 2 : 1,
        true,
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::lineCellFace::lineCellFace
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSet(name, mesh, dict),
    start_(dict.lookup("start")),
    end_(dict.lookup("end")),
    storeBothProcFaces_(false)
{}


Foam::sampledSets::lineCellFace::lineCellFace
(
    const word& name,
    const polyMesh& mesh,
    const word& axis,
    const point& start,
    const point& end
)
:
    sampledSet(name, mesh, axis),
    start_(start),
    end_(end),
    storeBothProcFaces_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::lineCellFace::~lineCellFace()
{}


// ************************************************************************* //
