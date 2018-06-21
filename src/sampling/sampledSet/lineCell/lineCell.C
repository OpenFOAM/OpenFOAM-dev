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

#include "lineCell.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(lineCell, 0);
    addToRunTimeSelectionTable(sampledSet, lineCell, word);
}
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::sampledSets::lineCell::calcMidPointSample
(
    const polyMesh& mesh,
    const point& prevPt,
    const label prevFace,
    const label prevSegment,
    const scalar prevCurveDist,
    const point& nextPt,
    const label nextFace,
    const label nextSegment,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    if (prevSegment == nextSegment)
    {
        const point pt = (prevPt + nextPt)/2;
        const vector delta = nextPt - prevPt;

        const label prevOwner = mesh.faceOwner()[prevFace];
        const label prevNeighbour =
            prevFace < mesh.faceNeighbour().size()
          ? mesh.faceNeighbour()[prevFace]
          : -1;
        const label nextOwner = mesh.faceOwner()[nextFace];
        const label nextNeighbour =
            nextFace < mesh.faceNeighbour().size()
          ? mesh.faceNeighbour()[nextFace]
          : -2;

        label celli = -1;
        if (prevOwner == nextOwner || prevOwner == nextNeighbour)
        {
            celli = prevOwner;
        }
        else if (prevNeighbour == nextOwner || prevNeighbour == nextNeighbour)
        {
            celli = prevNeighbour;
        }
        else
        {
            FatalErrorInFunction
                << "Adjacent faces in the same segment do not share a cell. "
                << "This is a bug." << exit(FatalError);
        }

        samplingPts.append(pt);
        samplingCells.append(celli);
        samplingFaces.append(-1);
        samplingCurveDist.append(prevCurveDist + mag(delta)/2);
        samplingSegments.append(prevSegment);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::lineCell::calcSamples
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
    lineFace::calcSamples
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

    // Append all mid points to the set
    for (label facei = 1; facei < facePts.size(); ++ facei)
    {
        calcMidPointSample
        (
            mesh(),
            facePts[facei - 1],
            faceFaces[facei - 1],
            faceSegments[facei - 1],
            faceCurveDist[facei - 1],
            facePts[facei],
            faceFaces[facei],
            faceSegments[facei],
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingSegments,
            samplingCurveDist
        );
    }
}


void Foam::sampledSets::lineCell::genSamples()
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

Foam::sampledSets::lineCell::lineCell
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::lineCell::~lineCell()
{}


// ************************************************************************* //
