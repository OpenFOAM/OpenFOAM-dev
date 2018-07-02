/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "points.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(points, 0);
    addToRunTimeSelectionTable(sampledSet, points, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::points::calcSamplesUnordered
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    forAll(points_, i)
    {
        const point& pt = points_[i];
        const label celli = searchEngine().findCell(pt);

        if (celli != -1)
        {
            samplingPts.append(pt);
            samplingCells.append(celli);
            samplingFaces.append(-1);
            samplingSegments.append(samplingSegments.size());
            samplingCurveDist.append(scalar(i));
        }
    }
}

void Foam::sampledSets::points::calcSamplesOrdered
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    const label n = points_.size();

    label sampleSegmentI = 0;
    label sampleI = 0;
    scalar sampleDist = 0;

    while (sampleI < n)
    {
        const point pt = points_[sampleI];

        const label sampleCellI = searchEngine().findCell(pt);

        if (sampleCellI == -1)
        {
            if (++ sampleI < n)
            {
                sampleDist += mag(points_[sampleI] - points_[sampleI - 1]);
            }
        }
        else
        {
            passiveParticle sampleParticle(mesh(), pt, sampleCellI);

            do
            {
                samplingPts.append(sampleParticle.position());
                samplingCells.append(sampleParticle.cell());
                samplingFaces.append(-1);
                samplingSegments.append(sampleSegmentI);
                samplingCurveDist.append(sampleDist);

                if (++ sampleI < n)
                {
                    const vector s = points_[sampleI] - points_[sampleI - 1];
                    sampleDist += mag(s);
                    sampleParticle.track(s, 0);
                }
            }
            while (sampleI < n && !sampleParticle.onBoundaryFace());

            ++ sampleSegmentI;
        }
    }
}


void Foam::sampledSets::points::genSamples()
{
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    if (!ordered_)
    {
        calcSamplesUnordered
        (
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingSegments,
            samplingCurveDist
        );
    }
    else
    {
        calcSamplesOrdered
        (
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingSegments,
            samplingCurveDist
        );
    }

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::points::points
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    points_(dict.lookup("points")),
    ordered_(dict.lookup("ordered"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::points::~points()
{}


// ************************************************************************* //
