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

#include "midPointAndFaceSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(midPointAndFaceSet, 0);
    addToRunTimeSelectionTable(sampledSet, midPointAndFaceSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::midPointAndFaceSet::genSamples()
{
    // Generate midpoints and add to face points

    List<point> mpfSamplePoints(3*size());
    labelList mpfSampleCells(3*size());
    labelList mpfSampleFaces(3*size());
    labelList mpfSampleSegments(3*size());
    scalarList mpfSampleCurveDist(3*size());

    label mpfSamplei = 0;
    label samplei = 0;

    while (size() > 0)
    {
        // Add first face
        mpfSamplePoints[mpfSamplei] = operator[](samplei);
        mpfSampleCells[mpfSamplei] = cells_[samplei];
        mpfSampleFaces[mpfSamplei] = faces_[samplei];
        mpfSampleSegments[mpfSamplei] = segments_[samplei];
        mpfSampleCurveDist[mpfSamplei] = curveDist_[samplei];
        mpfSamplei++;

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
                mpfSamplePoints[mpfSamplei] = midPoint;
                mpfSampleCells[mpfSamplei] = cellm;
                mpfSampleFaces[mpfSamplei] = -1;
                mpfSampleSegments[mpfSamplei] = segments_[samplei];
                mpfSampleCurveDist[mpfSamplei] =
                    mag(mpfSamplePoints[mpfSamplei] - start());

                mpfSamplei++;
            }

            // Add second face
            mpfSamplePoints[mpfSamplei] = operator[](samplei+1);
            mpfSampleCells[mpfSamplei] = cells_[samplei+1];
            mpfSampleFaces[mpfSamplei] = faces_[samplei+1];
            mpfSampleSegments[mpfSamplei] = segments_[samplei+1];
            mpfSampleCurveDist[mpfSamplei] =
                mag(mpfSamplePoints[mpfSamplei] - start());

            mpfSamplei++;

            samplei++;
        }

        if (samplei == size() - 1)
        {
            break;
        }
        samplei++;
    }

    mpfSamplePoints.setSize(mpfSamplei);
    mpfSampleCells.setSize(mpfSamplei);
    mpfSampleFaces.setSize(mpfSamplei);
    mpfSampleSegments.setSize(mpfSamplei);
    mpfSampleCurveDist.setSize(mpfSamplei);

    setSamples
    (
        mpfSamplePoints,
        mpfSampleCells,
        mpfSampleFaces,
        mpfSampleSegments,
        mpfSampleCurveDist
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::midPointAndFaceSet::midPointAndFaceSet
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


Foam::midPointAndFaceSet::midPointAndFaceSet
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

Foam::midPointAndFaceSet::~midPointAndFaceSet()
{}


// ************************************************************************* //
