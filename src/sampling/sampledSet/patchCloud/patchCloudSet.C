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

#include "patchCloudSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "Time.H"
#include "meshTools.H"
// For 'nearInfo' helper class only
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchCloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchCloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchCloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    if (debug)
    {
        Info<< "patchCloudSet : sampling on patches :" << endl;
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
    treeBoundBox bb(point::max, point::min);
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const polyPatch& pp = mesh().boundaryMesh()[iter.key()];

        forAll(pp, i)
        {
            patchFaces[sz++] = pp.start()+i;
        }

        // Do not do reduction.
        const boundBox patchBb(pp.points(), pp.meshPoints(), false);

        bb.min() = min(bb.min(), patchBb.min());
        bb.max() = max(bb.max(), patchBb.max());
    }
    bb = bb.extend(1e-4);

    indexedOctree<treeDataFace> patchTree
    (
        treeDataFace    // all information needed to search faces
        (
            false,      // do not cache bb
            mesh(),
            patchFaces  // boundary faces only
        ),
        bb,             // overall search domain
        8,              // maxLevel
        10,             // leafsize
        3.0             // duplicity
    );

    // Force calculation of face-diagonal decomposition
    (void)mesh().tetBasePtIs();


    // All the info for nearest. Construct to miss
    List<mappedPatchBase::nearInfo> nearest(sampleCoords_.size());

    forAll(sampleCoords_, sampleI)
    {
        const point& sample = sampleCoords_[sampleI];

        pointIndexHit& nearInfo = nearest[sampleI].first();

        // Find the nearest locally
        if (patchFaces.size())
        {
            nearInfo = patchTree.findNearest(sample, sqr(searchDist_));
        }
        else
        {
            nearInfo.setMiss();
        }


        // Fill in the distance field and the processor field
        if (!nearInfo.hit())
        {
            nearest[sampleI].second().first() = Foam::sqr(great);
            nearest[sampleI].second().second() = Pstream::myProcNo();
        }
        else
        {
            // Set nearest to mesh face label
            nearInfo.setIndex(patchFaces[nearInfo.index()]);

            nearest[sampleI].second().first() = magSqr
            (
                nearInfo.hitPoint()
              - sample
            );
            nearest[sampleI].second().second() = Pstream::myProcNo();
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);


    if (debug && Pstream::master())
    {
        OFstream str
        (
            mesh().time().path()
          / name()
          + "_nearest.obj"
        );
        Info<< "Dumping mapping as lines from supplied points to"
            << " nearest patch face to file " << str.name() << endl;

        label vertI = 0;

        forAll(nearest, i)
        {
            if (nearest[i].first().hit())
            {
                meshTools::writeOBJ(str, sampleCoords_[i]);
                vertI++;
                meshTools::writeOBJ(str, nearest[i].first().hitPoint());
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Store the sampling locations on the nearest processor
    forAll(nearest, sampleI)
    {
        const pointIndexHit& nearInfo = nearest[sampleI].first();

        if (nearInfo.hit())
        {
            if (nearest[sampleI].second().second() == Pstream::myProcNo())
            {
                label facei = nearInfo.index();

                samplingPts.append(nearInfo.hitPoint());
                samplingCells.append(mesh().faceOwner()[facei]);
                samplingFaces.append(facei);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0 * sampleI);
            }
        }
        else
        {
            // No processor found point near enough. Mark with special value
            // which is intercepted when interpolating
            if (Pstream::master())
            {
                samplingPts.append(sampleCoords_[sampleI]);
                samplingCells.append(-1);
                samplingFaces.append(-1);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0 * sampleI);
            }
        }
    }
}


void Foam::patchCloudSet::genSamples()
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

Foam::patchCloudSet::patchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords,
    const labelHashSet& patchSet,
    const scalar searchDist
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords),
    patchSet_(patchSet),
    searchDist_(searchDist)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::patchCloudSet::patchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.lookup("points")),
    patchSet_
    (
        mesh.boundaryMesh().patchSet
        (
            wordReList(dict.lookup("patches"))
        )
    ),
    searchDist_(readScalar(dict.lookup("maxDistance")))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchCloudSet::~patchCloudSet()
{}


// ************************************************************************* //
