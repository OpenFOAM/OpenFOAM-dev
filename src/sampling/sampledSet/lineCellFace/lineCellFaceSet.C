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

#include "lineCellFaceSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lineCellFaceSet, 0);
    addToRunTimeSelectionTable(sampledSet, lineCellFaceSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lineCellFaceSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Run the algorithm from lineFaceSet to get all the face intersections
    DynamicList<point> facePts;
    DynamicList<label> faceCells;
    DynamicList<label> faceFaces;
    DynamicList<label> faceSegments;
    DynamicList<scalar> faceCurveDist;
    lineFaceSet::calcSamples
    (
        mesh(),
        searchEngine(),
        start_,
        end_,
        facePts,
        faceCells,
        faceFaces,
        faceSegments,
        faceCurveDist
    );

    // If there are no intersections then quit
    if (!facePts.size())
    {
        return;
    }

    // Append all the face intersections to the set, additionally adding mid
    // points when the segment is the same
    samplingPts.append(facePts[0]);
    samplingCells.append(faceCells[0]);
    samplingFaces.append(faceFaces[0]);
    samplingSegments.append(faceSegments[0]);
    samplingCurveDist.append(faceCurveDist[0]);

    for (label facei = 1; facei < facePts.size(); ++ facei)
    {
        lineCellSet::calcMidPointSample
        (
            mesh(),
            samplingPts.last(),
            samplingFaces.last(),
            samplingSegments.last(),
            samplingCurveDist.last(),
            facePts[facei],
            faceFaces[facei],
            faceSegments[facei],
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingSegments,
            samplingCurveDist
        );

        samplingPts.append(facePts[facei]);
        samplingCells.append(faceCells[facei]);
        samplingFaces.append(faceFaces[facei]);
        samplingSegments.append(faceSegments[facei]);
        samplingCurveDist.append(faceCurveDist[facei]);
    }
}


void Foam::lineCellFaceSet::genSamples()
{
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

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

Foam::lineCellFaceSet::lineCellFaceSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    start_(dict.lookup("start")),
    end_(dict.lookup("end"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::lineCellFaceSet::lineCellFaceSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& start,
    const point& end
)
:
    sampledSet(name, mesh, searchEngine, axis),
    start_(start),
    end_(end)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lineCellFaceSet::~lineCellFaceSet()
{}


// ************************************************************************* //
