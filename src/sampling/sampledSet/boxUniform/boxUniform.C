/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "boxUniform.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(boxUniform, 0);
    addToRunTimeSelectionTable(sampledSet, boxUniform, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::boxUniform::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    for (label k = 0; k < nPoints_.z(); ++ k)
    {
        for (label j = 0; j < nPoints_.y(); ++ j)
        {
            for (label i = 0; i < nPoints_.x(); ++ i)
            {
                const vector t =
                    cmptDivide(vector(i, j, k), vector(nPoints_) - vector::one);

                const point pt =
                    cmptMultiply(vector::one - t, box_.min())
                  + cmptMultiply(t, box_.max());

                const label celli = searchEngine().findCell(pt);

                if (celli != -1)
                {
                    samplingPositions.append(pt);
                    samplingSegments.append
                    (
                        i + j*nPoints_.x() + k*nPoints_.x()*nPoints_.y()
                    );
                    samplingCells.append(celli);
                    samplingFaces.append(-1);
                }
            }
        }
    }
}


void Foam::sampledSets::boxUniform::genSamples()
{
    DynamicList<point> samplingPositions;
    DynamicList<label> samplingSegments;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;

    calcSamples
    (
        samplingPositions,
        samplingSegments,
        samplingCells,
        samplingFaces
    );

    samplingPositions.shrink();
    samplingSegments.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();

    setSamples
    (
        samplingPositions,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::boxUniform::boxUniform
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    box_(dict.lookup("box")),
    nPoints_(dict.lookup("nPoints"))
{
    genSamples();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::boxUniform::~boxUniform()
{}


// ************************************************************************* //
