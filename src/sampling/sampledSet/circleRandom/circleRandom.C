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

#include "circleRandom.H"
#include "meshSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(circleRandom, 0);
    addToRunTimeSelectionTable(sampledSet, circleRandom, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::circleRandom::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>&,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    const meshSearch& searchEngine = meshSearch::New(mesh());

    randomGenerator rndGen(261782, true);

    const vector radial1 = normalised(perpendicular(normal_));
    const vector radial2 = normalised(normal_ ^ radial1);

    for (label i = 0; i < nPoints_; ++ i)
    {
        // Request all random numbers simultaneously on all processors so that
        // the generator state stays consistent

        const scalar r = radius_*rndGen.scalar01();
        const scalar theta = 2*constant::mathematical::pi*rndGen.scalar01();
        const scalar c = cos(theta), s = sin(theta);

        const point pt = centre_ + r*(c*radial1 + s*radial2);
        const label celli = searchEngine.findCell(pt);

        if (celli != -1)
        {
            samplingPositions.append(pt);
            samplingSegments.append(i);
            samplingCells.append(celli);
            samplingFaces.append(-1);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::circleRandom::circleRandom
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSet(name, mesh, dict),
    centre_(dict.lookup("centre")),
    normal_(normalised(dict.lookup<vector>("normal"))),
    radius_(dict.lookup<scalar>("radius")),
    nPoints_(dict.lookup<label>("nPoints"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::circleRandom::~circleRandom()
{}


// ************************************************************************* //
