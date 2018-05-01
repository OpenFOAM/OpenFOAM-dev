/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "patchSeedSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "Time.H"
#include "meshTools.H"
//#include "Random.H"
// For 'facePoint' helper function only
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchSeedSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchSeedSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchSeedSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    if (debug)
    {
        Info<< "patchSeedSet : sampling on patches :" << endl;
    }

    // Construct search tree for all patch faces.
    label sz = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];

        sz += pp.size();

        if (debug)
        {
            Info<< "    " << pp.name() << " size " << pp.size() << endl;
        }
    }

    labelList patchFaces(sz);
    sz = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];
        forAll(pp, i)
        {
            patchFaces[sz++] = pp.start()+i;
        }
    }


    label totalSize = returnReduce(sz, sumOp<label>());


    // Shuffle and truncate if in random mode
    if (maxPoints_ < totalSize)
    {
        // Check what fraction of maxPoints_ I need to generate locally.
        label myMaxPoints = label(scalar(sz)/totalSize*maxPoints_);

        rndGenPtr_.reset(new Random(123456));
        Random& rndGen = rndGenPtr_();

        labelList subset = identity(sz);
        for (label iter = 0; iter < 4; iter++)
        {
            forAll(subset, i)
            {
                label j = rndGen.integer(0, subset.size()-1);
                Swap(subset[i], subset[j]);
            }
        }
        // Truncate
        subset.setSize(myMaxPoints);

        // Subset patchFaces
        patchFaces = UIndirectList<label>(patchFaces, subset)();

        if (debug)
        {
            Pout<< "In random mode : selected " << patchFaces.size()
                << " faces out of " << sz << endl;
        }
    }


    // Get points on patchFaces.
    globalIndex globalSampleNumbers(patchFaces.size());

    samplingPts.setCapacity(patchFaces.size());
    samplingCells.setCapacity(patchFaces.size());
    samplingFaces.setCapacity(patchFaces.size());
    samplingSegments.setCapacity(patchFaces.size());
    samplingCurveDist.setCapacity(patchFaces.size());

    // For calculation of min-decomp tet base points
    (void)mesh().tetBasePtIs();

    forAll(patchFaces, i)
    {
        label facei = patchFaces[i];
        pointIndexHit info = mappedPatchBase::facePoint
        (
            mesh(),
            facei,
            polyMesh::FACE_DIAG_TRIS
        );
        label celli = mesh().faceOwner()[facei];

        if (info.hit())
        {
            // Move the point into the cell
            const point& cc = mesh().cellCentres()[celli];
            samplingPts.append
            (
                info.hitPoint() + 1e-1*(cc-info.hitPoint())
            );
        }
        else
        {
            samplingPts.append(info.rawPoint());
        }
        samplingCells.append(celli);
        samplingFaces.append(facei);
        samplingSegments.append(0);
        samplingCurveDist.append(globalSampleNumbers.toGlobal(i));
    }
}


void Foam::patchSeedSet::genSamples()
{
    // Storage for sample points
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

Foam::patchSeedSet::patchSeedSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    patchSet_
    (
        mesh.boundaryMesh().patchSet
        (
            wordReList(dict.lookup("patches"))
        )
    ),
    // searchDist_(readScalar(dict.lookup("maxDistance"))),
    // offsetDist_(readScalar(dict.lookup("offsetDist"))),
    maxPoints_(readLabel(dict.lookup("maxPoints")))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchSeedSet::~patchSeedSet()
{}


// ************************************************************************* //
