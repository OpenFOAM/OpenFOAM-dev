/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "refinementFeatures.H"
#include "Time.H"
#include "Tuple2.H"
#include "DynamicField.H"
#include "featureEdgeMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementFeatures::read
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
{
    forAll(featDicts, featI)
    {
        const dictionary& dict = featDicts[featI];

        fileName featFileName(dict.lookup("file"));


        // Try reading extendedEdgeMesh first

        IOobject extFeatObj
        (
            featFileName,                       // name
            io.time().constant(),               // instance
            "extendedFeatureEdgeMesh",          // local
            io.time(),                          // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        const fileName fName(typeFilePath<extendedFeatureEdgeMesh>(extFeatObj));

        if (!fName.empty() && extendedEdgeMesh::canRead(fName))
        {
            autoPtr<extendedEdgeMesh> eMeshPtr = extendedEdgeMesh::New
            (
                fName
            );

            Info<< "Read extendedFeatureEdgeMesh " << extFeatObj.name()
                << nl << incrIndent;
            eMeshPtr().writeStats(Info);
            Info<< decrIndent << endl;

            set(featI, new extendedFeatureEdgeMesh(extFeatObj, eMeshPtr()));
        }
        else
        {
            // Try reading edgeMesh

            IOobject featObj
            (
                featFileName,
                io.time().constant(),
                searchableSurface::geometryDir(io.time()),
                io.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName fName(typeFilePath<featureEdgeMesh>(featObj));

            if (fName.empty())
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Could not open " << featObj.objectPath()
                    << exit(FatalIOError);
            }


            // Read as edgeMesh
            autoPtr<edgeMesh> eMeshPtr = edgeMesh::New(fName);
            const edgeMesh& eMesh = eMeshPtr();

            Info<< "Read edgeMesh " << featObj.name() << nl
                << incrIndent;
            eMesh.writeStats(Info);
            Info<< decrIndent << endl;


            // Analyse for feature points. These are all classified as mixed
            // points for lack of anything better
            const labelListList& pointEdges = eMesh.pointEdges();

            labelList oldToNew(eMesh.points().size(), -1);
            DynamicField<point> newPoints(eMesh.points().size());
            forAll(pointEdges, pointi)
            {
                if (pointEdges[pointi].size() > 2)
                {
                    oldToNew[pointi] = newPoints.size();
                    newPoints.append(eMesh.points()[pointi]);
                }
                // else if (pointEdges[pointi].size() == 2)
                // MEJ: do something based on a feature angle?
            }
            label nFeatures = newPoints.size();
            forAll(oldToNew, pointi)
            {
                if (oldToNew[pointi] == -1)
                {
                    oldToNew[pointi] = newPoints.size();
                    newPoints.append(eMesh.points()[pointi]);
                }
            }


            const edgeList& edges = eMesh.edges();
            edgeList newEdges(edges.size());
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                newEdges[edgeI] = edge
                (
                    oldToNew[e[0]],
                    oldToNew[e[1]]
                );
            }

            // Construct an extendedEdgeMesh with
            // - all points on more than 2 edges : mixed feature points
            // - all edges : external edges

            extendedEdgeMesh eeMesh
            (
                newPoints,          // pts
                newEdges,           // eds
                0,                  // (point) concaveStart
                0,                  // (point) mixedStart
                nFeatures,          // (point) nonFeatureStart
                edges.size(),       // (edge) internalStart
                edges.size(),       // (edge) flatStart
                edges.size(),       // (edge) openStart
                edges.size(),       // (edge) multipleStart
                vectorField(0),     // normals
                List<extendedEdgeMesh::sideVolumeType>(0),// normalVolumeTypes
                vectorField(0),     // edgeDirections
                labelListList(0),   // normalDirections
                labelListList(0),   // edgeNormals
                labelListList(0),   // featurePointNormals
                labelListList(0),   // featurePointEdges
                identity(newEdges.size())   // regionEdges
            );

            // Info<< "Constructed extendedFeatureEdgeMesh " << featObj.name()
            //    << nl << incrIndent;
            // eeMesh.writeStats(Info);
            // Info<< decrIndent << endl;

            set(featI, new extendedFeatureEdgeMesh(featObj, eeMesh));
        }

        const extendedEdgeMesh& eMesh = operator[](featI);

        if (dict.found("levels"))
        {
            List<Tuple2<scalar, label>> distLevels(dict["levels"]);

            if (dict.size() < 1)
            {
                FatalErrorInFunction
                    << " : levels should be at least size 1" << endl
                    << "levels : "  << dict["levels"]
                    << exit(FatalError);
            }

            distances_[featI].setSize(distLevels.size());
            levels_[featI].setSize(distLevels.size());

            forAll(distLevels, j)
            {
                distances_[featI][j] = distLevels[j].first();
                levels_[featI][j] = distLevels[j].second();

                // Check in incremental order
                if (j > 0)
                {
                    if
                    (
                        (distances_[featI][j] <= distances_[featI][j-1])
                     || (levels_[featI][j] > levels_[featI][j-1])
                    )
                    {
                        FatalErrorInFunction
                            << " : Refinement should be specified in order"
                            << " of increasing distance"
                            << " (and decreasing refinement level)." << endl
                            << "Distance:" << distances_[featI][j]
                            << " refinementLevel:" << levels_[featI][j]
                            << exit(FatalError);
                    }
                }
            }
        }
        else
        {
            // Look up 'level' for single level
            levels_[featI] = labelList(1, dict.lookup<label>("level"));
            distances_[featI] = scalarField(1, 0.0);
        }

        Info<< "Refinement level according to distance to "
            << featFileName << " (" << eMesh.points().size() << " points, "
            << eMesh.edges().size() << " edges)." << endl;
        forAll(levels_[featI], j)
        {
            Info<< "    level " << levels_[featI][j]
                << " for all cells within " << distances_[featI][j]
                << " metre." << endl;
        }
    }
}


void Foam::refinementFeatures::buildTrees(const label featI)
{
    const extendedEdgeMesh& eMesh = operator[](featI);
    const pointField& points = eMesh.points();
    const edgeList& edges = eMesh.edges();

    // Calculate bb of all points
    treeBoundBox bb(points);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    bb = bb.extend(1e-4);

    edgeTrees_.set
    (
        featI,
        new indexedOctree<treeDataEdge>
        (
            treeDataEdge
            (
                false,                  // do not cache bb
                edges,
                points,
                identity(edges.size())
            ),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );


    labelList featurePoints(identity(eMesh.nonFeatureStart()));

    pointTrees_.set
    (
        featI,
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(points, featurePoints),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );
}


// Find maximum level of a shell.
void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const label featI,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[featI];

    const scalarField& distances = distances_[featI];

    // Collect all those points that have a current maxLevel less than
    // (any of) the shell. Also collect the furthest distance allowable
    // to any shell with a higher level.

    pointField candidates(pt.size());
    labelList candidateMap(pt.size());
    scalarField candidateDistSqr(pt.size());
    label candidateI = 0;

    forAll(maxLevel, pointi)
    {
        forAllReverse(levels, levelI)
        {
            if (levels[levelI] > maxLevel[pointi])
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateDistSqr[candidateI] = sqr(distances[levelI]);
                candidateI++;
                break;
            }
        }
    }
    candidates.setSize(candidateI);
    candidateMap.setSize(candidateI);
    candidateDistSqr.setSize(candidateI);

    // Do the expensive nearest test only for the candidate points.
    const indexedOctree<treeDataEdge>& tree = edgeTrees_[featI];

    List<pointIndexHit> nearInfo(candidates.size());
    forAll(candidates, candidateI)
    {
        nearInfo[candidateI] = tree.findNearest
        (
            candidates[candidateI],
            candidateDistSqr[candidateI]
        );
    }

    // Update maxLevel
    forAll(nearInfo, candidateI)
    {
        if (nearInfo[candidateI].hit())
        {
            // Check which level it actually is in.
            label minDistI = findLower
            (
                distances,
                mag(nearInfo[candidateI].hitPoint()-candidates[candidateI])
            );

            label pointi = candidateMap[candidateI];

            // pt is in between shell[minDistI] and shell[minDistI+1]
            maxLevel[pointi] = levels[minDistI+1];
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge>>&
Foam::refinementFeatures::regionEdgeTrees() const
{
    if (!regionEdgeTreesPtr_.valid())
    {
        regionEdgeTreesPtr_.reset
        (
            new PtrList<indexedOctree<treeDataEdge>>(size())
        );
        PtrList<indexedOctree<treeDataEdge>>& trees = regionEdgeTreesPtr_();

        forAll(*this, featI)
        {
            const extendedEdgeMesh& eMesh = operator[](featI);
            const pointField& points = eMesh.points();
            const edgeList& edges = eMesh.edges();

            // Calculate bb of all points
            treeBoundBox bb(points);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(1e-4);

            trees.set
            (
                featI,
                new indexedOctree<treeDataEdge>
                (
                    treeDataEdge
                    (
                        false,                  // do not cache bb
                        edges,
                        points,
                        eMesh.regionEdges()
                    ),
                    bb,     // overall search domain
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                )
            );
        }
    }
    return regionEdgeTreesPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementFeatures::refinementFeatures
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
:
    PtrList<extendedFeatureEdgeMesh>(featDicts.size()),
    distances_(featDicts.size()),
    levels_(featDicts.size()),
    edgeTrees_(featDicts.size()),
    pointTrees_(featDicts.size())
{
    // Read features
    read(io, featDicts);

    // Search engines
    forAll(*this, i)
    {
        buildTrees(i);
    }
}


//Foam::refinementFeatures::refinementFeatures
//(
//    const objectRegistry& io,
//    const PtrList<dictionary>& featDicts,
//    const scalar minCos
//)
//:
//    PtrList<extendedFeatureEdgeMesh>(featDicts.size()),
//    distances_(featDicts.size()),
//    levels_(featDicts.size()),
//    edgeTrees_(featDicts.size()),
//    pointTrees_(featDicts.size())
//{
//    // Read features
//    read(io, featDicts);
//
//    // Search engines
//    forAll(*this, i)
//    {
//        const edgeMesh& eMesh = operator[](i);
//        const pointField& points = eMesh.points();
//        const edgeList& edges = eMesh.edges();
//        const labelListList& pointEdges = eMesh.pointEdges();
//
//        DynamicList<label> featurePoints;
//        forAll(pointEdges, pointi)
//        {
//            const labelList& pEdges = pointEdges[pointi];
//            if (pEdges.size() > 2)
//            {
//                featurePoints.append(pointi);
//            }
//            else if (pEdges.size() == 2)
//            {
//                // Check the angle
//                const edge& e0 = edges[pEdges[0]];
//                const edge& e1 = edges[pEdges[1]];
//
//                const point& p = points[pointi];
//                const point& p0 = points[e0.otherVertex(pointi)];
//                const point& p1 = points[e1.otherVertex(pointi)];
//
//                vector v0 = p-p0;
//                scalar v0Mag = mag(v0);
//
//                vector v1 = p1-p;
//                scalar v1Mag = mag(v1);
//
//                if
//                (
//                    v0Mag > small
//                 && v1Mag > small
//                 && ((v0/v0Mag & v1/v1Mag) < minCos)
//                )
//                {
//                    featurePoints.append(pointi);
//                }
//            }
//        }
//
//        Info<< "Detected " << featurePoints.size()
//            << " featurePoints out of " << points.size()
//            << " points on feature " << i   // eMesh.name()
//            << " when using feature cos " << minCos << endl;
//
//        buildTrees(i, featurePoints);
//    }
//}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::refinementFeatures::findNearestEdge
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo,
    vectorField& nearNormal
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());
    nearInfo = pointIndexHit();
    nearNormal.setSize(samples.size());
    nearNormal = Zero;

    forAll(edgeTrees_, featI)
    {
        const indexedOctree<treeDataEdge>& tree = edgeTrees_[featI];

        if (tree.shapes().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearInfo[sampleI].hit())
                {
                    distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[sampleI] = featI;
                    nearInfo[sampleI] = pointIndexHit
                    (
                        info.hit(),
                        info.hitPoint(),
                        tree.shapes().edgeLabels()[info.index()]
                    );

                    const treeDataEdge& td = tree.shapes();
                    const edge& e = td.edges()[nearInfo[sampleI].index()];
                    nearNormal[sampleI] =  e.vec(td.points());
                    nearNormal[sampleI] /= mag(nearNormal[sampleI])+vSmall;
                }
            }
        }
    }
}


void Foam::refinementFeatures::findNearestRegionEdge
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo,
    vectorField& nearNormal
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());
    nearInfo = pointIndexHit();
    nearNormal.setSize(samples.size());
    nearNormal = Zero;


    const PtrList<indexedOctree<treeDataEdge>>& regionTrees =
        regionEdgeTrees();

    forAll(regionTrees, featI)
    {
        const indexedOctree<treeDataEdge>& regionTree = regionTrees[featI];

        forAll(samples, sampleI)
        {
            const point& sample = samples[sampleI];

            scalar distSqr;
            if (nearInfo[sampleI].hit())
            {
                distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
            }
            else
            {
                distSqr = nearestDistSqr[sampleI];
            }

            // Find anything closer than current best
            pointIndexHit info = regionTree.findNearest(sample, distSqr);

            if (info.hit())
            {
                const treeDataEdge& td = regionTree.shapes();

                nearFeature[sampleI] = featI;
                nearInfo[sampleI] = pointIndexHit
                (
                    info.hit(),
                    info.hitPoint(),
                    regionTree.shapes().edgeLabels()[info.index()]
                );

                const edge& e = td.edges()[nearInfo[sampleI].index()];
                nearNormal[sampleI] =  e.vec(td.points());
                nearNormal[sampleI] /= mag(nearNormal[sampleI])+vSmall;
            }
        }
    }
}


//void Foam::refinementFeatures::findNearestPoint
//(
//    const pointField& samples,
//    const scalarField& nearestDistSqr,
//    labelList& nearFeature,
//    labelList& nearIndex
//) const
//{
//    nearFeature.setSize(samples.size());
//    nearFeature = -1;
//    nearIndex.setSize(samples.size());
//    nearIndex = -1;
//
//    forAll(pointTrees_, featI)
//    {
//        const indexedOctree<treeDataPoint>& tree = pointTrees_[featI];
//
//        if (tree.shapes().pointLabels().size() > 0)
//        {
//            forAll(samples, sampleI)
//            {
//                const point& sample = samples[sampleI];
//
//                scalar distSqr;
//                if (nearFeature[sampleI] != -1)
//                {
//                    label nearFeatI = nearFeature[sampleI];
//                    const indexedOctree<treeDataPoint>& nearTree =
//                        pointTrees_[nearFeatI];
//                    label featPointi =
//                        nearTree.shapes().pointLabels()[nearIndex[sampleI]];
//                    const point& featPt =
//                        operator[](nearFeatI).points()[featPointi];
//                    distSqr = magSqr(featPt-sample);
//                }
//                else
//                {
//                    distSqr = nearestDistSqr[sampleI];
//                }
//
//                pointIndexHit info = tree.findNearest(sample, distSqr);
//
//                if (info.hit())
//                {
//                    nearFeature[sampleI] = featI;
//                    nearIndex[sampleI] = info.index();
//                }
//            }
//        }
//    }
//}


void Foam::refinementFeatures::findNearestPoint
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());
    nearInfo = pointIndexHit();

    forAll(pointTrees_, featI)
    {
        const indexedOctree<treeDataPoint>& tree = pointTrees_[featI];

        if (tree.shapes().pointLabels().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearFeature[sampleI] != -1)
                {
                    distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[sampleI] = featI;
                    nearInfo[sampleI] = pointIndexHit
                    (
                        info.hit(),
                        info.hitPoint(),
                        tree.shapes().pointLabels()[info.index()]
                    );
                }
            }
        }
    }
}


void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(*this, featI)
    {
        findHigherLevel(pt, featI, maxLevel);
    }
}


Foam::scalar Foam::refinementFeatures::maxDistance() const
{
    scalar overallMax = -great;
    forAll(distances_, featI)
    {
        overallMax = max(overallMax, max(distances_[featI]));
    }
    return overallMax;
}


// ************************************************************************* //
