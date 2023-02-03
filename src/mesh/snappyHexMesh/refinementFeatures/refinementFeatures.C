/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    forAll(featDicts, feati)
    {
        const dictionary& dict = featDicts[feati];

        fileName featFileName(dict.lookup("file"));


        // Try reading extendedEdgeMesh first

        typeIOobject<extendedFeatureEdgeMesh> extFeatObj
        (
            featFileName,                       // name
            io.time().constant(),               // instance
            "extendedFeatureEdgeMesh",          // local
            io.time(),                          // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        const fileName fName(extFeatObj.filePath());

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

            set(feati, new extendedFeatureEdgeMesh(extFeatObj, eMeshPtr()));
        }
        else
        {
            // Try reading edgeMesh

            typeIOobject<featureEdgeMesh> featObj
            (
                featFileName,
                io.time().constant(),
                searchableSurface::geometryDir(io.time()),
                io.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName fName(featObj.filePath());

            if (fName.empty())
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Could not open "
                    << featObj.objectPath()
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
            forAll(edges, edgei)
            {
                const edge& e = edges[edgei];
                newEdges[edgei] = edge
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
                identityMap(newEdges.size())   // regionEdges
            );

            // Info<< "Constructed extendedFeatureEdgeMesh " << featObj.name()
            //    << nl << incrIndent;
            // eeMesh.writeStats(Info);
            // Info<< decrIndent << endl;

            set(feati, new extendedFeatureEdgeMesh(featObj, eeMesh));
        }

        const extendedEdgeMesh& eMesh = operator[](feati);

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

            distances_[feati].setSize(distLevels.size());
            levels_[feati].setSize(distLevels.size());

            forAll(distLevels, j)
            {
                distances_[feati][j] = distLevels[j].first();
                levels_[feati][j] = distLevels[j].second();

                // Check in incremental order
                if (j > 0)
                {
                    if
                    (
                        (distances_[feati][j] <= distances_[feati][j-1])
                     || (levels_[feati][j] > levels_[feati][j-1])
                    )
                    {
                        FatalErrorInFunction
                            << " : Refinement should be specified in order"
                            << " of increasing distance"
                            << " (and decreasing refinement level)." << endl
                            << "Distance:" << distances_[feati][j]
                            << " refinementLevel:" << levels_[feati][j]
                            << exit(FatalError);
                    }
                }
            }
        }
        else
        {
            // Look up 'level' for single level
            levels_[feati] = labelList(1, dict.lookup<label>("level"));
            distances_[feati] = scalarField(1, 0.0);
        }

        Info<< "Refinement level according to distance to "
            << featFileName << " (" << eMesh.points().size() << " points, "
            << eMesh.edges().size() << " edges)." << endl;
        forAll(levels_[feati], j)
        {
            Info<< "    level " << levels_[feati][j]
                << " for all cells within " << distances_[feati][j]
                << " metre." << endl;
        }
    }
}


void Foam::refinementFeatures::buildTrees(const label feati)
{
    const extendedEdgeMesh& eMesh = operator[](feati);
    const pointField& points = eMesh.points();
    const edgeList& edges = eMesh.edges();

    // Calculate bb of all points
    treeBoundBox bb(points);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    bb = bb.extend(1e-4);

    edgeTrees_.set
    (
        feati,
        new indexedOctree<treeDataEdge>
        (
            treeDataEdge
            (
                false,                  // do not cache bb
                edges,
                points,
                identityMap(edges.size())
            ),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        )
    );


    labelList featurePoints(identityMap(eMesh.nonFeatureStart()));

    pointTrees_.set
    (
        feati,
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


// Find maximum level of a feature edge.
void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const label feati,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[feati];

    const scalarField& distances = distances_[feati];

    // Collect all those points that have a current maxLevel less than
    // (any of) the feature edge. Also collect the furthest distance allowable
    // to any feature edge with a higher level.

    pointField candidates(pt.size());
    labelList candidateMap(pt.size());
    scalarField candidateDistSqr(pt.size());
    label candidatei = 0;

    forAll(maxLevel, pointi)
    {
        forAllReverse(levels, leveli)
        {
            if (levels[leveli] > maxLevel[pointi])
            {
                candidates[candidatei] = pt[pointi];
                candidateMap[candidatei] = pointi;
                candidateDistSqr[candidatei] = sqr(distances[leveli]);
                candidatei++;
                break;
            }
        }
    }
    candidates.setSize(candidatei);
    candidateMap.setSize(candidatei);
    candidateDistSqr.setSize(candidatei);

    // Do the expensive nearest test only for the candidate points.
    const indexedOctree<treeDataEdge>& tree = edgeTrees_[feati];

    List<pointIndexHit> nearInfo(candidates.size());
    forAll(candidates, candidatei)
    {
        nearInfo[candidatei] = tree.findNearest
        (
            candidates[candidatei],
            candidateDistSqr[candidatei]
        );
    }

    // Update maxLevel
    forAll(nearInfo, candidatei)
    {
        if (nearInfo[candidatei].hit())
        {
            // Check which level it actually is in.
            label minDistI = findLower
            (
                distances,
                mag(nearInfo[candidatei].hitPoint()-candidates[candidatei])
            );

            label pointi = candidateMap[candidatei];

            // pt is in between feature[minDistI] and feature[minDistI+1]
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

        forAll(*this, feati)
        {
            const extendedEdgeMesh& eMesh = operator[](feati);
            const pointField& points = eMesh.points();
            const edgeList& edges = eMesh.edges();

            // Calculate bb of all points
            treeBoundBox bb(points);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(1e-4);

            trees.set
            (
                feati,
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

    forAll(edgeTrees_, feati)
    {
        const indexedOctree<treeDataEdge>& tree = edgeTrees_[feati];

        if (tree.shapes().size() > 0)
        {
            forAll(samples, samplei)
            {
                const point& sample = samples[samplei];

                scalar distSqr;
                if (nearInfo[samplei].hit())
                {
                    distSqr = magSqr(nearInfo[samplei].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[samplei];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[samplei] = feati;
                    nearInfo[samplei] = pointIndexHit
                    (
                        info.hit(),
                        info.hitPoint(),
                        tree.shapes().edgeLabels()[info.index()]
                    );

                    const treeDataEdge& td = tree.shapes();
                    const edge& e = td.edges()[nearInfo[samplei].index()];
                    nearNormal[samplei] =  e.vec(td.points());
                    nearNormal[samplei] /= mag(nearNormal[samplei])+vSmall;
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

    forAll(regionTrees, feati)
    {
        const indexedOctree<treeDataEdge>& regionTree = regionTrees[feati];

        forAll(samples, samplei)
        {
            const point& sample = samples[samplei];

            scalar distSqr;
            if (nearInfo[samplei].hit())
            {
                distSqr = magSqr(nearInfo[samplei].hitPoint()-sample);
            }
            else
            {
                distSqr = nearestDistSqr[samplei];
            }

            // Find anything closer than current best
            pointIndexHit info = regionTree.findNearest(sample, distSqr);

            if (info.hit())
            {
                const treeDataEdge& td = regionTree.shapes();

                nearFeature[samplei] = feati;
                nearInfo[samplei] = pointIndexHit
                (
                    info.hit(),
                    info.hitPoint(),
                    regionTree.shapes().edgeLabels()[info.index()]
                );

                const edge& e = td.edges()[nearInfo[samplei].index()];
                nearNormal[samplei] =  e.vec(td.points());
                nearNormal[samplei] /= mag(nearNormal[samplei])+vSmall;
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
//    forAll(pointTrees_, feati)
//    {
//        const indexedOctree<treeDataPoint>& tree = pointTrees_[feati];
//
//        if (tree.shapes().pointLabels().size() > 0)
//        {
//            forAll(samples, samplei)
//            {
//                const point& sample = samples[samplei];
//
//                scalar distSqr;
//                if (nearFeature[samplei] != -1)
//                {
//                    label nearFeatI = nearFeature[samplei];
//                    const indexedOctree<treeDataPoint>& nearTree =
//                        pointTrees_[nearFeatI];
//                    label featPointi =
//                        nearTree.shapes().pointLabels()[nearIndex[samplei]];
//                    const point& featPt =
//                        operator[](nearFeatI).points()[featPointi];
//                    distSqr = magSqr(featPt-sample);
//                }
//                else
//                {
//                    distSqr = nearestDistSqr[samplei];
//                }
//
//                pointIndexHit info = tree.findNearest(sample, distSqr);
//
//                if (info.hit())
//                {
//                    nearFeature[samplei] = feati;
//                    nearIndex[samplei] = info.index();
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

    forAll(pointTrees_, feati)
    {
        const indexedOctree<treeDataPoint>& tree = pointTrees_[feati];

        if (tree.shapes().pointLabels().size() > 0)
        {
            forAll(samples, samplei)
            {
                const point& sample = samples[samplei];

                scalar distSqr;
                if (nearFeature[samplei] != -1)
                {
                    distSqr = magSqr(nearInfo[samplei].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[samplei];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[samplei] = feati;
                    nearInfo[samplei] = pointIndexHit
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
    // Maximum level of any feature edge. Start off with level of point.
    maxLevel = ptLevel;

    forAll(*this, feati)
    {
        findHigherLevel(pt, feati, maxLevel);
    }
}


Foam::scalar Foam::refinementFeatures::maxDistance() const
{
    scalar overallMax = -great;
    forAll(distances_, feati)
    {
        overallMax = max(overallMax, max(distances_[feati]));
    }
    return overallMax;
}


// ************************************************************************* //
