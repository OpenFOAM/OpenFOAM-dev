/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "midPointSet.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(midPointSet, 0);
    addToRunTimeSelectionTable(sampledSet, midPointSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::midPointSet::genSamples()
{
    // Generate midpoints.

    List<point> midPoints(2*size());
    labelList midCells(2*size());
    labelList midSegments(2*size());
    scalarList midCurveDist(2*size());

    label mSamplei = 0;
    label samplei = 0;

    while (size() > 0)
    {
        // Calculate midpoint between samplei and samplei+1 (if in same segment)
        while
        (
            (samplei < size() - 1)
         && (segments_[samplei] == segments_[samplei+1])
        )
        {
            point midPoint(0.5*(operator[](samplei) + operator[](samplei+1)));
            label cellm = pointInCell(midPoint, samplei);

            if (cellm != -1)
            {
                midPoints[mSamplei] = midPoint;
                midCells[mSamplei] = cellm;
                midSegments[mSamplei] = segments_[samplei];
                midCurveDist[mSamplei] = mag(midPoints[mSamplei] - start());
                mSamplei++;
            }

            samplei++;
        }

        if (samplei == size() - 1)
        {
            break;
        }

        samplei++;
    }

    midPoints.setSize(mSamplei);
    midCells.setSize(mSamplei);
    midSegments.setSize(mSamplei);
    midCurveDist.setSize(mSamplei);

    setSamples
    (
        midPoints,
        midCells,
        labelList(midCells.size(), -1),
        midSegments,
        midCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::midPointSet::midPointSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& start,
    const point& end
)
:
    faceOnlySet(name, mesh, searchEngine, axis, start, end)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::midPointSet::midPointSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    faceOnlySet(name, mesh, searchEngine, dict)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::midPointSet::~midPointSet()
{}


// ************************************************************************* //
