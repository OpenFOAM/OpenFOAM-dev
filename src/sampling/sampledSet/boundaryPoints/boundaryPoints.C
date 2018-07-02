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

#include "boundaryPoints.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "Time.H"
#include "meshTools.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(boundaryPoints, 0);
    addToRunTimeSelectionTable(sampledSet, boundaryPoints, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::boundaryPoints::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Construct a single list of all patch faces
    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patches_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];
        nPatchFaces += pp.size();
    }
    labelList patchFaces(nPatchFaces);
    nPatchFaces = 0;
    forAllConstIter(labelHashSet, patches_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];
        forAll(pp, i)
        {
            patchFaces[nPatchFaces++] = pp.start()+i;
        }
    }

    // Construct a processor-local bound box
    treeBoundBox patchBB(point::max, point::min);
    forAllConstIter(labelHashSet, patches_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];
        const boundBox patchBb(pp.points(), pp.meshPoints(), false);
        patchBB.min() = min(patchBB.min(), patchBb.min());
        patchBB.max() = max(patchBB.max(), patchBb.max());
    }
    patchBB = patchBB.extend(1e-4);

    // Create a search tree
    indexedOctree<treeDataFace> patchTree
    (
        treeDataFace    // all information needed to search faces
        (
            false,      // do not cache the bound box
            mesh(),
            patchFaces
        ),
        patchBB,        // overall search box
        8,              // maximum number of levels
        10,             // how many elements per leaf
        3.0             // in how many leaves is a shape on average
    );

    // Force calculation of face-diagonal decomposition
    (void)mesh().tetBasePtIs();

    // Generate the nearest patch information for each sampling point
    List<mappedPatchBase::nearInfo> nearest(points_.size());
    forAll(points_, sampleI)
    {
        const point& sample = points_[sampleI];

        pointIndexHit& nearHit = nearest[sampleI].first();
        scalar& nearDist = nearest[sampleI].second().first();
        label& nearProc = nearest[sampleI].second().second();

        // Find the nearest
        if (patchFaces.size())
        {
            nearHit = patchTree.findNearest(sample, sqr(maxDistance_));
        }
        else
        {
            nearHit.setMiss();
        }

        // Fill in the information
        if (nearHit.hit())
        {
            nearHit.setIndex(patchFaces[nearHit.index()]);
            nearDist = magSqr(nearHit.hitPoint() - sample);
            nearProc = Pstream::myProcNo();
        }
        else
        {
            nearHit.setIndex(-1);
            nearDist = Foam::sqr(great);
            nearProc = Pstream::myProcNo();
        }
    }

    // Reduce to get the nearest patch locations globally
    Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);

    // Dump connecting lines from the sampling points to the hit locations
    if (debug && Pstream::master())
    {
        OFstream str(mesh().time().path() / name() + "_nearest.obj");

        label verti = 0;

        forAll(nearest, i)
        {
            if (nearest[i].first().hit())
            {
                meshTools::writeOBJ(str, points_[i]);
                verti++;
                meshTools::writeOBJ(str, nearest[i].first().hitPoint());
                verti++;
                str << "l " << verti - 1 << ' ' << verti << nl;
            }
        }
    }

    // Store the sampling locations on the nearest processor
    forAll(nearest, sampleI)
    {
        const pointIndexHit& nearHit = nearest[sampleI].first();
        const label& nearProc = nearest[sampleI].second().second();

        if (nearHit.hit())
        {
            if (nearProc == Pstream::myProcNo())
            {
                label facei = nearHit.index();

                samplingPts.append(nearHit.hitPoint());
                samplingCells.append(mesh().faceOwner()[facei]);
                samplingFaces.append(facei);
                samplingSegments.append(0);
                samplingCurveDist.append(sampleI);
            }
        }
        else
        {
            WarningInFunction
                << "Unable to find location on patches " << patches_
                << " for the point " << points_[sampleI]
                << " within a distance of " << maxDistance_ << endl;
        }
    }
}


void Foam::sampledSets::boundaryPoints::genSamples()
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

Foam::sampledSets::boundaryPoints::boundaryPoints
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    points_(dict.lookup("points")),
    patches_
    (
        mesh.boundaryMesh().patchSet
        (
            wordReList(dict.lookup("patches"))
        )
    ),
    maxDistance_(readScalar(dict.lookup("maxDistance")))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::boundaryPoints::~boundaryPoints()
{}


// ************************************************************************* //
