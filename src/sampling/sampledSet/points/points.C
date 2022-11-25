/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "sampledSetCloud.H"
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


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::sampledSets::points::calcSamples
(
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const pointField& points,
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>& samplingDistances,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
)
{
    // Create a cloud with which to track segments
    sampledSetCloud particles
    (
        mesh,
        points::typeName,
        IDLList<sampledSetParticle>()
    );

    // Consider each point
    label segmenti = 0, samplei = 0, pointi0 = labelMax, pointi = 0;
    scalar distance = 0;
    while (pointi < points.size())
    {
        // Sum the distance to the start of the track
        for (label pointj = pointi0; pointj < pointi; ++ pointj)
        {
            distance += mag(points[pointj + 1] - points[pointj]);
        }

        // Update the old point index
        pointi0 = pointi;

        // Get unique processor and cell that this sample point is in
        const labelPair procAndCelli = returnReduce
        (
            labelPair
            (
                Pstream::myProcNo(),
                searchEngine.findCell(points[pointi])
            ),
            [](const labelPair& a, const labelPair& b)
            {
                return
                    a.second() != -1 && b.second() != -1
                  ? a.first() < b.first() ? a : b
                  : a.second() != -1 ? a : b;
            }
        );

        // Skip this point if it is not in the global mesh
        if (procAndCelli.second() == -1)
        {
            ++ pointi;
        }

        // If the point is in the global mesh then track to create a segment
        else
        {
            sampledSetParticle::trackingData td
            (
                particles,
                points,
                true,
                false,
                false,
                samplingPositions,
                samplingDistances,
                samplingCells,
                samplingFaces
            );

            // Clear the cloud, then, if the point is in this local mesh,
            // initialise a particle at the point
            particles.clear();
            if (procAndCelli.first() == Pstream::myProcNo())
            {
                particles.addParticle
                (
                    new sampledSetParticle
                    (
                        mesh,
                        points[pointi],
                        procAndCelli.second(),
                        pointi,
                        1,
                        distance
                    )
                );

                particles.first()->store(particles, td);
            }

            // Track to create this segment
            particles.move(particles, td);

            // Set the segment indices
            samplingSegments.append
            (
                labelList
                (
                    samplingPositions.size() - samplingSegments.size(),
                    segmenti
                )
            );

            // Move on to the next segment
            ++ segmenti;

            // Determine the global number of samples completed
            const label samplei0 = samplei;
            samplei = returnReduce(samplingPositions.size(), sumOp<label>());

            // Move to the next unsampled point
            pointi += samplei - samplei0;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::points::calcSamplesUnordered
(
    DynamicList<point>& samplingPositions,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    forAll(points_, i)
    {
        const point& pt = points_[i];
        const label celli = searchEngine().findCell(pt);

        if (celli != -1)
        {
            samplingPositions.append(pt);
            samplingSegments.append(i);
            samplingCells.append(celli);
            samplingFaces.append(-1);
        }
    }
}


void Foam::sampledSets::points::calcSamplesOrdered
(
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>& samplingDistances,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    // Calculate the sampling topology
    calcSamples
    (
        mesh(),
        searchEngine(),
        pointField(points_),
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


void Foam::sampledSets::points::genSamples()
{
    DynamicList<point> samplingPositions;
    DynamicList<scalar> samplingDistances;
    DynamicList<label> samplingSegments;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;

    if (!ordered_)
    {
        calcSamplesUnordered
        (
            samplingPositions,
            samplingSegments,
            samplingCells,
            samplingFaces
        );
    }
    else
    {
        calcSamplesOrdered
        (
            samplingPositions,
            samplingDistances,
            samplingSegments,
            samplingCells,
            samplingFaces
        );
    }

    samplingPositions.shrink();
    samplingDistances.shrink();
    samplingSegments.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();

    if (!ordered_)
    {
        setSamples
        (
            samplingPositions,
            samplingSegments,
            samplingCells,
            samplingFaces
        );
    }
    else
    {
        setSamples
        (
            samplingPositions,
            samplingDistances,
            samplingSegments,
            samplingCells,
            samplingFaces
        );
    }
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::points::~points()
{}


// ************************************************************************* //
