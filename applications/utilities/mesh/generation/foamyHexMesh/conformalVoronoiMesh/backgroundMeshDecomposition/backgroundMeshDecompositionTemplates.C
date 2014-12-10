/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "backgroundMeshDecomposition.H"
#include "pointConversion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<typename PointType>
Foam::autoPtr<Foam::mapDistribute>
Foam::backgroundMeshDecomposition::distributePoints
(
    List<PointType>& points
) const
{
    labelList toProc(processorPosition(points));

    autoPtr<mapDistribute> map(buildMap(toProc));

    map().distribute(points);

    return map;
}


template<typename PointType>
Foam::labelList Foam::backgroundMeshDecomposition::processorPosition
(
    const List<PointType>& pts
) const
{
    DynamicList<label> toCandidateProc;
    DynamicList<point> testPoints;
    labelList ptBlockStart(pts.size(), -1);
    labelList ptBlockSize(pts.size(), -1);

    label nTotalCandidates = 0;

    forAll(pts, pI)
    {
        const pointFromPoint pt = topoint(pts[pI]);

        label nCandidates = 0;

        forAll(allBackgroundMeshBounds_, procI)
        {
            if (allBackgroundMeshBounds_[procI].contains(pt))
            {
                toCandidateProc.append(procI);
                testPoints.append(pt);

                nCandidates++;
            }
        }

        ptBlockStart[pI] = nTotalCandidates;
        ptBlockSize[pI] = nCandidates;

        nTotalCandidates += nCandidates;
    }

    // Needed for reverseDistribute
    label preDistributionToCandidateProcSize = toCandidateProc.size();

    autoPtr<mapDistribute> map(buildMap(toCandidateProc));

    map().distribute(testPoints);

    List<bool> pointOnCandidate(testPoints.size(), false);

    // Test candidate points on candidate processors
    forAll(testPoints, tPI)
    {
        pointOnCandidate[tPI] = positionOnThisProcessor(testPoints[tPI]);
    }

    map().reverseDistribute
    (
        preDistributionToCandidateProcSize,
        pointOnCandidate
    );

    labelList ptProc(pts.size(), -1);

    DynamicList<label> failedPointIndices;
    DynamicList<point> failedPoints;

    forAll(pts, pI)
    {
        // Extract the sub list of results for this point

        SubList<bool> ptProcResults
        (
            pointOnCandidate,
            ptBlockSize[pI],
            ptBlockStart[pI]
        );

        forAll(ptProcResults, pPRI)
        {
            if (ptProcResults[pPRI])
            {
                ptProc[pI] = toCandidateProc[ptBlockStart[pI] + pPRI];

                break;
            }
        }

        if (ptProc[pI] < 0)
        {
            const pointFromPoint pt = topoint(pts[pI]);

            if (!globalBackgroundBounds_.contains(pt))
            {
                FatalErrorIn
                (
                    "Foam::labelList"
                    "Foam::backgroundMeshDecomposition::processorPosition"
                    "("
                        "const List<point>&"
                    ") const"
                )
                    << "The position " << pt
                    << " is not in any part of the background mesh "
                    << globalBackgroundBounds_ << endl
                    << "A background mesh with a wider margin around "
                    << "the geometry may help."
                    << exit(FatalError);
            }

            if (debug)
            {
                WarningIn
                (
                    "Foam::labelList"
                    "Foam::backgroundMeshDecomposition::processorPosition"
                    "("
                        "const List<point>&"
                    ") const"
                )   << "The position " << pt
                    << " was not found in the background mesh "
                    << globalBackgroundBounds_ << ", finding nearest."
                    << endl;
            }

            failedPointIndices.append(pI);
            failedPoints.append(pt);
        }
    }

    labelList ptNearestProc(processorNearestPosition(failedPoints));

    forAll(failedPoints, fPI)
    {
        label pI = failedPointIndices[fPI];

        ptProc[pI] = ptNearestProc[fPI];
    }

    return ptProc;
}


// ************************************************************************* //
