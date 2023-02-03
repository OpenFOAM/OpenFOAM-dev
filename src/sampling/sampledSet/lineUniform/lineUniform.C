/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "lineUniform.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "sampledSetCloud.H"
#include "points.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(lineUniform, 0);
    addToRunTimeSelectionTable(sampledSet, lineUniform, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::lineUniform::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>& samplingDistances,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    // Calculate all sampling points
    const scalarField ts(scalarList(identityMap(nPoints_))/(nPoints_ - 1));
    const pointField points((1 - ts)*start_ + ts*end_);

    // Calculate the sampling topology
    points::calcSamples
    (
        mesh(),
        searchEngine(),
        points,
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );

    // Overwrite the distances
    forAll(samplingPositions, i)
    {
        samplingDistances[i] = mag(samplingPositions[i] - start_);
    }
}


void Foam::sampledSets::lineUniform::genSamples()
{
    DynamicList<point> samplingPositions;
    DynamicList<scalar> samplingDistances;
    DynamicList<label> samplingSegments;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;

    calcSamples
    (
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );

    samplingPositions.shrink();
    samplingDistances.shrink();
    samplingSegments.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();

    setSamples
    (
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::lineUniform::lineUniform
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    start_(dict.lookup("start")),
    end_(dict.lookup("end")),
    nPoints_(dict.lookup<label>("nPoints"))
{
    genSamples();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::lineUniform::~lineUniform()
{}


// ************************************************************************* //
