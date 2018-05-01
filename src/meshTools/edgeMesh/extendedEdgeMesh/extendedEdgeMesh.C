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

#include "extendedEdgeMesh.H"
#include "surfaceFeatures.H"
#include "triSurface.H"
#include "Random.H"
#include "Time.H"
#include "OBJstream.H"
#include "DynamicField.H"
#include "edgeMeshFormatsCore.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedEdgeMesh, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::extendedEdgeMesh::pointStatus,
        4
    >::names[] =
    {
        "convex",
        "concave",
        "mixed",
        "nonFeature"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::extendedEdgeMesh::edgeStatus,
        6
    >::names[] =
    {
        "external",
        "internal",
        "flat",
        "open",
        "multiple",
        "none"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::extendedEdgeMesh::sideVolumeType,
        4
    >::names[] =
    {
        "inside",
        "outside",
        "both",
        "neither"
    };
}

const Foam::NamedEnum<Foam::extendedEdgeMesh::pointStatus, 4>
    Foam::extendedEdgeMesh::pointStatusNames_;

const Foam::NamedEnum<Foam::extendedEdgeMesh::edgeStatus, 6>
    Foam::extendedEdgeMesh::edgeStatusNames_;

const Foam::NamedEnum<Foam::extendedEdgeMesh::sideVolumeType, 4>
    Foam::extendedEdgeMesh::sideVolumeTypeNames_;

Foam::scalar Foam::extendedEdgeMesh::cosNormalAngleTol_ =
    Foam::cos(degToRad(0.1));


Foam::label Foam::extendedEdgeMesh::convexStart_ = 0;


Foam::label Foam::extendedEdgeMesh::externalStart_ = 0;


Foam::label Foam::extendedEdgeMesh::nPointTypes = 4;


Foam::label Foam::extendedEdgeMesh::nEdgeTypes = 5;


Foam::wordHashSet Foam::extendedEdgeMesh::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


Foam::wordHashSet Foam::extendedEdgeMesh::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}



// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::extendedEdgeMesh::canReadType
(
    const word& ext,
    const bool verbose
)
{
    return edgeMeshFormatsCore::checkSupport
    (
        readTypes(),
        ext,
        verbose,
        "reading"
   );
}


bool Foam::extendedEdgeMesh::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return edgeMeshFormatsCore::checkSupport
    (
        writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


bool Foam::extendedEdgeMesh::canRead
(
    const fileName& name,
    const bool verbose
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::extendedEdgeMesh::pointStatus
Foam::extendedEdgeMesh::classifyFeaturePoint
(
    label ptI
) const
{
    const labelList& ptEds(pointEdges()[ptI]);

    label nPtEds = ptEds.size();
    label nExternal = 0;
    label nInternal = 0;

    if (nPtEds == 0)
    {
        // There are no edges attached to the point, this is a problem
        return NONFEATURE;
    }

    forAll(ptEds, i)
    {
        edgeStatus edStat = getEdgeStatus(ptEds[i]);

        if (edStat == EXTERNAL)
        {
            nExternal++;
        }
        else if (edStat == INTERNAL)
        {
            nInternal++;
        }
    }

    if (nExternal == nPtEds)
    {
        return CONVEX;
    }
    else if (nInternal == nPtEds)
    {
        return CONCAVE;
    }
    else
    {
        return MIXED;
    }
}


Foam::extendedEdgeMesh::edgeStatus
Foam::extendedEdgeMesh::classifyEdge
(
    const List<vector>& norms,
    const labelList& edNorms,
    const vector& fC0tofC1
)
{
    label nEdNorms = edNorms.size();

    if (nEdNorms == 1)
    {
        return OPEN;
    }
    else if (nEdNorms == 2)
    {
        const vector& n0(norms[edNorms[0]]);
        const vector& n1(norms[edNorms[1]]);

        if ((n0 & n1) > cosNormalAngleTol_)
        {
            return FLAT;
        }
        else if ((fC0tofC1 & n0) > 0.0)
        {
            return INTERNAL;
        }
        else
        {
            return EXTERNAL;
        }
    }
    else if (nEdNorms > 2)
    {
        return MULTIPLE;
    }
    else
    {
        // There is a problem - the edge has no normals
        return NONE;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedEdgeMesh::extendedEdgeMesh()
:
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh(const extendedEdgeMesh& fem)
:
    edgeMesh(fem),
    concaveStart_(fem.concaveStart()),
    mixedStart_(fem.mixedStart()),
    nonFeatureStart_(fem.nonFeatureStart()),
    internalStart_(fem.internalStart()),
    flatStart_(fem.flatStart()),
    openStart_(fem.openStart()),
    multipleStart_(fem.multipleStart()),
    normals_(fem.normals()),
    normalVolumeTypes_(fem.normalVolumeTypes()),
    edgeDirections_(fem.edgeDirections()),
    normalDirections_(fem.normalDirections()),
    edgeNormals_(fem.edgeNormals()),
    featurePointNormals_(fem.featurePointNormals()),
    featurePointEdges_(fem.featurePointEdges()),
    regionEdges_(fem.regionEdges()),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh(Istream& is)
{
    is >> *this;
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
:
    edgeMesh(pointLst, edgeLst),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const surfaceFeatures& sFeat,
    const boolList& surfBaffleRegions
)
:
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(-1),
    mixedStart_(-1),
    nonFeatureStart_(-1),
    internalStart_(-1),
    flatStart_(-1),
    openStart_(-1),
    multipleStart_(-1),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    // Extract and reorder the data from surfaceFeatures
    const triSurface& surf = sFeat.surface();
    const labelList& featureEdges = sFeat.featureEdges();
    const labelList& featurePoints = sFeat.featurePoints();

    // Get a labelList of all the featureEdges that are region edges
    const labelList regionFeatureEdges(identity(sFeat.nRegionEdges()));

    sortPointsAndEdges
    (
        surf,
        featureEdges,
        regionFeatureEdges,
        featurePoints
    );

    const labelListList& edgeFaces = surf.edgeFaces();

    normalVolumeTypes_.setSize(normals_.size());

    // Noting when the normal of a face has been used so not to duplicate
    labelList faceMap(surf.size(), -1);

    label nAdded = 0;

    forAll(featureEdges, i)
    {
        label sFEI = featureEdges[i];

        // Pick up the faces adjacent to the feature edge
        const labelList& eFaces = edgeFaces[sFEI];

        forAll(eFaces, j)
        {
            label eFI = eFaces[j];

            // Check to see if the points have been already used
            if (faceMap[eFI] == -1)
            {
                normalVolumeTypes_[nAdded++] =
                    (
                        surfBaffleRegions[surf[eFI].region()]
                      ? BOTH
                      : INSIDE
                    );

                faceMap[eFI] = nAdded - 1;
            }
        }
    }
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const PrimitivePatch<face, List, pointField, point>& surf,
    const labelList& featureEdges,
    const labelList& regionFeatureEdges,
    const labelList& featurePoints
)
:
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(-1),
    mixedStart_(-1),
    nonFeatureStart_(-1),
    internalStart_(-1),
    flatStart_(-1),
    openStart_(-1),
    multipleStart_(-1),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    sortPointsAndEdges
    (
        surf,
        featureEdges,
        regionFeatureEdges,
        featurePoints
    );
}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const pointField& pts,
    const edgeList& eds,
    label concaveStart,
    label mixedStart,
    label nonFeatureStart,
    label internalStart,
    label flatStart,
    label openStart,
    label multipleStart,
    const vectorField& normals,
    const List<sideVolumeType>& normalVolumeTypes,
    const vectorField& edgeDirections,
    const labelListList& normalDirections,
    const labelListList& edgeNormals,
    const labelListList& featurePointNormals,
    const labelListList& featurePointEdges,
    const labelList& regionEdges
)
:
    edgeMesh(pts, eds),
    concaveStart_(concaveStart),
    mixedStart_(mixedStart),
    nonFeatureStart_(nonFeatureStart),
    internalStart_(internalStart),
    flatStart_(flatStart),
    openStart_(openStart),
    multipleStart_(multipleStart),
    normals_(normals),
    normalVolumeTypes_(normalVolumeTypes),
    edgeDirections_(edgeDirections),
    normalDirections_(normalDirections),
    edgeNormals_(edgeNormals),
    featurePointNormals_(featurePointNormals),
    featurePointEdges_(featurePointEdges),
    regionEdges_(regionEdges),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{}


Foam::extendedEdgeMesh::extendedEdgeMesh
(
    const fileName& name,
    const word& ext
)
:
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    read(name, ext);
}


Foam::extendedEdgeMesh::extendedEdgeMesh(const fileName& name)
:
    edgeMesh(pointField(0), edgeList(0)),
    concaveStart_(0),
    mixedStart_(0),
    nonFeatureStart_(0),
    internalStart_(0),
    flatStart_(0),
    openStart_(0),
    multipleStart_(0),
    normals_(0),
    normalVolumeTypes_(0),
    edgeDirections_(0),
    normalDirections_(0),
    edgeNormals_(0),
    featurePointNormals_(0),
    featurePointEdges_(0),
    regionEdges_(0),
    pointTree_(),
    edgeTree_(),
    edgeTreesByType_()
{
    read(name);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedEdgeMesh::~extendedEdgeMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::extendedEdgeMesh::read(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }
    else
    {
        return read(name, ext);
    }
}


// Read from file in given format
bool Foam::extendedEdgeMesh::read
(
    const fileName& name,
    const word& ext
)
{
    // read via selector mechanism
    transfer(New(name, ext)());
    return true;
}


void Foam::extendedEdgeMesh::nearestFeaturePoint
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = pointTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedEdgeMesh::nearestFeatureEdge
(
    const point& sample,
    scalar searchDistSqr,
    pointIndexHit& info
) const
{
    info = edgeTree().findNearest
    (
        sample,
        searchDistSqr
    );
}


void Foam::extendedEdgeMesh::nearestFeatureEdge
(
    const pointField& samples,
    const scalarField& searchDistSqr,
    pointIndexHitList& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        nearestFeatureEdge
        (
            samples[i],
            searchDistSqr[i],
            info[i]
        );
    }
}


void Foam::extendedEdgeMesh::nearestFeatureEdgeByType
(
    const point& sample,
    const scalarField& searchDistSqr,
    pointIndexHitList& info
) const
{
    const PtrList<indexedOctree<treeDataEdge>>& edgeTrees = edgeTreesByType();

    info.setSize(edgeTrees.size());

    labelList sliceStarts(edgeTrees.size());

    sliceStarts[0] = externalStart_;
    sliceStarts[1] = internalStart_;
    sliceStarts[2] = flatStart_;
    sliceStarts[3] = openStart_;
    sliceStarts[4] = multipleStart_;

    forAll(edgeTrees, i)
    {
        info[i] = edgeTrees[i].findNearest
        (
            sample,
            searchDistSqr[i]
        );

        // The index returned by the indexedOctree is local to the slice of
        // edges it was supplied with, return the index to the value in the
        // complete edge list

        info[i].setIndex(info[i].index() + sliceStarts[i]);
    }
}


void Foam::extendedEdgeMesh::allNearestFeaturePoints
(
    const point& sample,
    scalar searchRadiusSqr,
    pointIndexHitList& info
) const
{
    // Pick up all the feature points that intersect the search sphere
    labelList elems = pointTree().findSphere
    (
        sample,
        searchRadiusSqr
    );

    DynamicList<pointIndexHit> dynPointHit(elems.size());

    forAll(elems, elemI)
    {
        label index = elems[elemI];
        label ptI = pointTree().shapes().pointLabels()[index];
        const point& pt = points()[ptI];

        pointIndexHit nearHit(true, pt, index);

        dynPointHit.append(nearHit);
    }

    info.transfer(dynPointHit);
}


void Foam::extendedEdgeMesh::allNearestFeatureEdges
(
    const point& sample,
    const scalar searchRadiusSqr,
    pointIndexHitList& info
) const
{
    const PtrList<indexedOctree<treeDataEdge>>& edgeTrees = edgeTreesByType();

    info.setSize(edgeTrees.size());

    labelList sliceStarts(edgeTrees.size());

    sliceStarts[0] = externalStart_;
    sliceStarts[1] = internalStart_;
    sliceStarts[2] = flatStart_;
    sliceStarts[3] = openStart_;
    sliceStarts[4] = multipleStart_;

    DynamicList<pointIndexHit> dynEdgeHit(edgeTrees.size()*3);

    // Loop over all the feature edge types
    forAll(edgeTrees, i)
    {
        // Pick up all the edges that intersect the search sphere
        labelList elems = edgeTrees[i].findSphere
        (
            sample,
            searchRadiusSqr
        );

        forAll(elems, elemI)
        {
            label index = elems[elemI];
            label edgeI = edgeTrees[i].shapes().edgeLabels()[index];
            const edge& e = edges()[edgeI];

            pointHit hitPoint = e.line(points()).nearestDist(sample);

            label hitIndex = index + sliceStarts[i];

            pointIndexHit nearHit
            (
                hitPoint.hit(),
                hitPoint.rawPoint(),
                hitIndex
            );

            dynEdgeHit.append(nearHit);
        }
    }

    info.transfer(dynEdgeHit);
}


Foam::scalar Foam::extendedEdgeMesh::minDisconnectedDist
(
    const pointIndexHitList& hitList
) const
{
    scalar minDist = GREAT;

    for
    (
        label hi1 = 0;
        hi1 < hitList.size() - 1;
        ++hi1
    )
    {
        const pointIndexHit& pHit1 = hitList[hi1];

        if (pHit1.hit())
        {
            const edge& e1 = edges()[pHit1.index()];

            for
            (
                label hi2 = hi1 + 1;
                hi2 < hitList.size();
                ++hi2
            )
            {
                const pointIndexHit& pHit2 = hitList[hi2];

                if (pHit2.hit())
                {
                    const edge& e2 = edges()[pHit2.index()];

                    // Don't refine if the edges are connected to each other
                    if (!e1.connected(e2))
                    {
                        scalar curDist =
                            mag(pHit1.hitPoint() - pHit2.hitPoint());

                        minDist = min(curDist, minDist);
                    }
                }
            }
        }
    }

    return minDist;
}


const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::extendedEdgeMesh::pointTree() const
{
    if (pointTree_.empty())
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4)
        );

        bb.min() -= point(rootVSmall, rootVSmall, rootVSmall);
        bb.max() += point(rootVSmall, rootVSmall, rootVSmall);

        const labelList featurePointLabels = identity(nonFeatureStart_);

        pointTree_.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint
                (
                    points(),
                    featurePointLabels
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return pointTree_();
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::extendedEdgeMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4)
        );

        bb.min() -= point(rootVSmall, rootVSmall, rootVSmall);
        bb.max() += point(rootVSmall, rootVSmall, rootVSmall);

        labelList allEdges(identity(edges().size()));

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    points(),       // points
                    allEdges        // selected edges
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return edgeTree_();
}


const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge>>&
Foam::extendedEdgeMesh::edgeTreesByType() const
{
    if (edgeTreesByType_.size() == 0)
    {
        edgeTreesByType_.setSize(nEdgeTypes);

        Random rndGen(872141);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1e-4)
        );

        bb.min() -= point(rootVSmall, rootVSmall, rootVSmall);
        bb.max() += point(rootVSmall, rootVSmall, rootVSmall);

        labelListList sliceEdges(nEdgeTypes);

        // External edges
        sliceEdges[0] =
            identity(internalStart_ - externalStart_) + externalStart_;

        // Internal edges
        sliceEdges[1] = identity(flatStart_ - internalStart_) + internalStart_;

        // Flat edges
        sliceEdges[2] = identity(openStart_ - flatStart_) + flatStart_;

        // Open edges
        sliceEdges[3] = identity(multipleStart_ - openStart_) + openStart_;

        // Multiple edges
        sliceEdges[4] =
            identity(edges().size() - multipleStart_) + multipleStart_;

        forAll(edgeTreesByType_, i)
        {
            edgeTreesByType_.set
            (
                i,
                new indexedOctree<treeDataEdge>
                (
                    treeDataEdge
                    (
                        false,          // cachebb
                        edges(),        // edges
                        points(),       // points
                        sliceEdges[i]   // selected edges
                    ),
                    bb,     // bb
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                )
            );
        }
    }

    return edgeTreesByType_;
}


void Foam::extendedEdgeMesh::transfer(extendedEdgeMesh& mesh)
{
    edgeMesh::transfer(mesh);

    concaveStart_ = mesh.concaveStart_;
    mixedStart_ = mesh.mixedStart_;
    nonFeatureStart_ = mesh.nonFeatureStart_;
    internalStart_ = mesh.internalStart_;
    flatStart_ = mesh.flatStart_;
    openStart_ = mesh.openStart_;
    multipleStart_ = mesh.multipleStart_;
    normals_.transfer(mesh.normals_);
    normalVolumeTypes_.transfer(mesh.normalVolumeTypes_);
    edgeDirections_.transfer(mesh.edgeDirections_);
    normalDirections_.transfer(mesh.normalDirections_);
    edgeNormals_.transfer(mesh.edgeNormals_);
    featurePointNormals_.transfer(mesh.featurePointNormals_);
    featurePointEdges_.transfer(mesh.featurePointEdges_);
    regionEdges_.transfer(mesh.regionEdges_);
    pointTree_ = mesh.pointTree_;
    edgeTree_ = mesh.edgeTree_;
    edgeTreesByType_.transfer(mesh.edgeTreesByType_);
}


Foam::Xfer<Foam::extendedEdgeMesh> Foam::extendedEdgeMesh::xfer()
{
    return xferMoveTo<extendedEdgeMesh, extendedEdgeMesh>(*this);
}


void Foam::extendedEdgeMesh::clear()
{
    edgeMesh::clear();
    concaveStart_ = 0;
    mixedStart_ = 0;
    nonFeatureStart_ = 0;
    internalStart_ = 0;
    flatStart_ = 0;
    openStart_ = 0;
    multipleStart_ = 0;
    normals_.clear();
    normalVolumeTypes_.clear();
    edgeDirections_.clear();
    normalDirections_.clear();
    edgeNormals_.clear();
    featurePointNormals_.clear();
    featurePointEdges_.clear();
    regionEdges_.clear();
    pointTree_.clear();
    edgeTree_.clear();
    edgeTreesByType_.clear();
}


void Foam::extendedEdgeMesh::add(const extendedEdgeMesh& fem)
{
    // Points
    // ~~~~~~

    // From current points into combined points
    labelList reversePointMap(points().size());
    labelList reverseFemPointMap(fem.points().size());

    label newPointi = 0;
    for (label i = 0; i < concaveStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = 0; i < fem.concaveStart(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    // Concave
    label newConcaveStart = newPointi;
    for (label i = concaveStart(); i < mixedStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = fem.concaveStart(); i < fem.mixedStart(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    // Mixed
    label newMixedStart = newPointi;
    for (label i = mixedStart(); i < nonFeatureStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = fem.mixedStart(); i < fem.nonFeatureStart(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    // Non-feature
    label newNonFeatureStart = newPointi;
    for (label i = nonFeatureStart(); i < points().size(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    for (label i = fem.nonFeatureStart(); i < fem.points().size(); i++)
    {
        reverseFemPointMap[i] = newPointi++;
    }

    pointField newPoints(newPointi);
    newPoints.rmap(points(), reversePointMap);
    newPoints.rmap(fem.points(), reverseFemPointMap);


    // Edges
    // ~~~~~

    // From current edges into combined edges
    labelList reverseEdgeMap(edges().size());
    labelList reverseFemEdgeMap(fem.edges().size());

    // External
    label newEdgeI = 0;
    for (label i = 0; i < internalStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = 0; i < fem.internalStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Internal
    label newInternalStart = newEdgeI;
    for (label i = internalStart(); i < flatStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.internalStart(); i < fem.flatStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Flat
    label newFlatStart = newEdgeI;
    for (label i = flatStart(); i < openStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.flatStart(); i < fem.openStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Open
    label newOpenStart = newEdgeI;
    for (label i = openStart(); i < multipleStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.openStart(); i < fem.multipleStart(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    // Multiple
    label newMultipleStart = newEdgeI;
    for (label i = multipleStart(); i < edges().size(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    for (label i = fem.multipleStart(); i < fem.edges().size(); i++)
    {
        reverseFemEdgeMap[i] = newEdgeI++;
    }

    edgeList newEdges(newEdgeI);
    forAll(edges(), i)
    {
        const edge& e = edges()[i];
        newEdges[reverseEdgeMap[i]] = edge
        (
            reversePointMap[e[0]],
            reversePointMap[e[1]]
        );
    }
    forAll(fem.edges(), i)
    {
        const edge& e = fem.edges()[i];
        newEdges[reverseFemEdgeMap[i]] = edge
        (
            reverseFemPointMap[e[0]],
            reverseFemPointMap[e[1]]
        );
    }

    pointField newEdgeDirections(newEdgeI);
    newEdgeDirections.rmap(edgeDirections(), reverseEdgeMap);
    newEdgeDirections.rmap(fem.edgeDirections(), reverseFemEdgeMap);




    // Normals
    // ~~~~~~~

    // Combine normals
    DynamicField<point> newNormals(normals().size()+fem.normals().size());
    newNormals.append(normals());
    newNormals.append(fem.normals());


    // Combine and re-index into newNormals
    labelListList newEdgeNormals(edgeNormals().size()+fem.edgeNormals().size());
    UIndirectList<labelList>(newEdgeNormals, reverseEdgeMap) =
        edgeNormals();
    UIndirectList<labelList>(newEdgeNormals, reverseFemEdgeMap) =
        fem.edgeNormals();
    forAll(reverseFemEdgeMap, i)
    {
        label mapI = reverseFemEdgeMap[i];
        labelList& en = newEdgeNormals[mapI];
        forAll(en, j)
        {
            en[j] += normals().size();
        }
    }


    // Combine and re-index into newFeaturePointNormals
    labelListList newFeaturePointNormals
    (
       featurePointNormals().size()
     + fem.featurePointNormals().size()
    );

    // Note: featurePointNormals only go up to nonFeatureStart
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reversePointMap, featurePointNormals().size())
    ) = featurePointNormals();
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reverseFemPointMap, fem.featurePointNormals().size())
    ) = fem.featurePointNormals();
    forAll(fem.featurePointNormals(), i)
    {
        label mapI = reverseFemPointMap[i];
        labelList& fn = newFeaturePointNormals[mapI];
        forAll(fn, j)
        {
            fn[j] += normals().size();
        }
    }


    // Combine regionEdges
    DynamicList<label> newRegionEdges
    (
        regionEdges().size()
      + fem.regionEdges().size()
    );
    forAll(regionEdges(), i)
    {
        newRegionEdges.append(reverseEdgeMap[regionEdges()[i]]);
    }
    forAll(fem.regionEdges(), i)
    {
        newRegionEdges.append(reverseFemEdgeMap[fem.regionEdges()[i]]);
    }


    // Assign
    // ~~~~~~

    // Transfer
    concaveStart_ = newConcaveStart;
    mixedStart_ = newMixedStart;
    nonFeatureStart_ = newNonFeatureStart;

    // Reset points and edges
    reset(xferMove(newPoints), newEdges.xfer());

    // Transfer
    internalStart_ = newInternalStart;
    flatStart_ = newFlatStart;
    openStart_ = newOpenStart;
    multipleStart_ = newMultipleStart;

    edgeDirections_.transfer(newEdgeDirections);

    normals_.transfer(newNormals);
    edgeNormals_.transfer(newEdgeNormals);
    featurePointNormals_.transfer(newFeaturePointNormals);

    regionEdges_.transfer(newRegionEdges);

    pointTree_.clear();
    edgeTree_.clear();
    edgeTreesByType_.clear();
}


void Foam::extendedEdgeMesh::flipNormals()
{
    // Points
    // ~~~~~~

    // From current points into new points
    labelList reversePointMap(identity(points().size()));

    // Flip convex and concave points

    label newPointi = 0;
    // Concave points become convex
    for (label i = concaveStart(); i < mixedStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }
    // Convex points become concave
    label newConcaveStart = newPointi;
    for (label i = 0; i < concaveStart(); i++)
    {
        reversePointMap[i] = newPointi++;
    }


    // Edges
    // ~~~~~~

    // From current edges into new edges
    labelList reverseEdgeMap(identity(edges().size()));

    // Flip external and internal edges

    label newEdgeI = 0;
    // Internal become external
    for (label i = internalStart(); i < flatStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }
    // External become internal
    label newInternalStart = newEdgeI;
    for (label i = 0; i < internalStart(); i++)
    {
        reverseEdgeMap[i] = newEdgeI++;
    }


    pointField newPoints(points().size());
    newPoints.rmap(points(), reversePointMap);

    edgeList newEdges(edges().size());
    forAll(edges(), i)
    {
        const edge& e = edges()[i];
        newEdges[reverseEdgeMap[i]] = edge
        (
            reversePointMap[e[0]],
            reversePointMap[e[1]]
        );
    }


    // Normals are flipped
    // ~~~~~~~~~~~~~~~~~~~

    pointField newEdgeDirections(edges().size());
    newEdgeDirections.rmap(-1.0*edgeDirections(), reverseEdgeMap);

    pointField newNormals(-1.0*normals());

    labelListList newEdgeNormals(edgeNormals().size());
    UIndirectList<labelList>(newEdgeNormals, reverseEdgeMap) = edgeNormals();

    labelListList newFeaturePointNormals(featurePointNormals().size());

    // Note: featurePointNormals only go up to nonFeatureStart
    UIndirectList<labelList>
    (
        newFeaturePointNormals,
        SubList<label>(reversePointMap, featurePointNormals().size())
    ) = featurePointNormals();

    labelList newRegionEdges(regionEdges().size());
    forAll(regionEdges(), i)
    {
        newRegionEdges[i] = reverseEdgeMap[regionEdges()[i]];
    }

    // Transfer
    concaveStart_ = newConcaveStart;

    // Reset points and edges
    reset(xferMove(newPoints), newEdges.xfer());

    // Transfer
    internalStart_ = newInternalStart;

    edgeDirections_.transfer(newEdgeDirections);
    normals_.transfer(newNormals);
    edgeNormals_.transfer(newEdgeNormals);
    featurePointNormals_.transfer(newFeaturePointNormals);
    regionEdges_.transfer(newRegionEdges);

    pointTree_.clear();
    edgeTree_.clear();
    edgeTreesByType_.clear();
}


void Foam::extendedEdgeMesh::writeObj
(
    const fileName& prefix,
    const bool verbose
) const
{
    Info<< nl << "Writing extendedEdgeMesh components to " << prefix
        << endl;

    edgeMesh::write(prefix + "_edgeMesh.obj");

    if (!verbose) return;

    OBJstream convexFtPtStr(prefix + "_convexFeaturePts.obj");
    Info<< "Writing convex feature points to " << convexFtPtStr.name() << endl;

    for(label i = 0; i < concaveStart_; i++)
    {
        convexFtPtStr.write(points()[i]);
    }

    OBJstream concaveFtPtStr(prefix + "_concaveFeaturePts.obj");
    Info<< "Writing concave feature points to "
        << concaveFtPtStr.name() << endl;

    for(label i = concaveStart_; i < mixedStart_; i++)
    {
        convexFtPtStr.write(points()[i]);
    }

    OBJstream mixedFtPtStr(prefix + "_mixedFeaturePts.obj");
    Info<< "Writing mixed feature points to " << mixedFtPtStr.name() << endl;

    for(label i = mixedStart_; i < nonFeatureStart_; i++)
    {
        mixedFtPtStr.write(points()[i]);
    }

    OBJstream mixedFtPtStructureStr(prefix + "_mixedFeaturePtsStructure.obj");
    Info<< "Writing mixed feature point structure to "
        << mixedFtPtStructureStr.name() << endl;

    for(label i = mixedStart_; i < nonFeatureStart_; i++)
    {
        const labelList& ptEds = pointEdges()[i];

        forAll(ptEds, j)
        {
            const edge& e = edges()[ptEds[j]];
            mixedFtPtStructureStr.write
            (
                linePointRef(points()[e[0]],
                points()[e[1]])
            );
        }
    }

    OBJstream externalStr(prefix + "_externalEdges.obj");
    Info<< "Writing external edges to " << externalStr.name() << endl;

    for (label i = externalStart_; i < internalStart_; i++)
    {
        const edge& e = edges()[i];
        externalStr.write(linePointRef(points()[e[0]], points()[e[1]]));
    }

    OBJstream internalStr(prefix + "_internalEdges.obj");
    Info<< "Writing internal edges to " << internalStr.name() << endl;

    for (label i = internalStart_; i < flatStart_; i++)
    {
        const edge& e = edges()[i];
        internalStr.write(linePointRef(points()[e[0]], points()[e[1]]));
    }

    OBJstream flatStr(prefix + "_flatEdges.obj");
    Info<< "Writing flat edges to " << flatStr.name() << endl;

    for (label i = flatStart_; i < openStart_; i++)
    {
        const edge& e = edges()[i];
        flatStr.write(linePointRef(points()[e[0]], points()[e[1]]));
    }

    OBJstream openStr(prefix + "_openEdges.obj");
    Info<< "Writing open edges to " << openStr.name() << endl;

    for (label i = openStart_; i < multipleStart_; i++)
    {
        const edge& e = edges()[i];
        openStr.write(linePointRef(points()[e[0]], points()[e[1]]));
    }

    OBJstream multipleStr(prefix + "_multipleEdges.obj");
    Info<< "Writing multiple edges to " << multipleStr.name() << endl;

    for (label i = multipleStart_; i < edges().size(); i++)
    {
        const edge& e = edges()[i];
        multipleStr.write(linePointRef(points()[e[0]], points()[e[1]]));
    }

    OBJstream regionStr(prefix + "_regionEdges.obj");
    Info<< "Writing region edges to " << regionStr.name() << endl;

    forAll(regionEdges_, i)
    {
        const edge& e = edges()[regionEdges_[i]];
        regionStr.write(linePointRef(points()[e[0]], points()[e[1]]));
    }

    OBJstream edgeDirsStr(prefix + "_edgeDirections.obj");
    Info<< "Writing edge directions to " << edgeDirsStr.name() << endl;

    forAll(edgeDirections_, i)
    {
        const vector& eVec = edgeDirections_[i];
        const edge& e = edges()[i];

        edgeDirsStr.write
        (
            linePointRef(points()[e.start()], eVec + points()[e.start()])
        );
    }
}


void Foam::extendedEdgeMesh::writeStats(Ostream& os) const
{
    edgeMesh::writeStats(os);

    os  << indent << "point classification :" << nl;
    os  << incrIndent;
    os  << indent << "convex feature points          : "
        << setw(8) << concaveStart_-convexStart_
        << nl;
    os  << indent << "concave feature points         : "
        << setw(8) << mixedStart_-concaveStart_
        << nl;
    os  << indent << "mixed feature points           : "
        << setw(8) << nonFeatureStart_-mixedStart_
        << nl;
    os  << indent << "other (non-feature) points     : "
        << setw(8) << points().size()-nonFeatureStart_
        << nl;
    os  << decrIndent;

    os  << indent << "edge classification :" << nl;
    os  << incrIndent;
    os  << indent << "external (convex angle) edges  : "
        << setw(8) << internalStart_-externalStart_
        << nl;
    os  << indent << "internal (concave angle) edges : "
        << setw(8) << flatStart_-internalStart_
        << nl;
    os  << indent << "flat region edges              : "
        << setw(8) << openStart_-flatStart_
        << nl;
    os  << indent << "open edges                     : "
        << setw(8) << multipleStart_-openStart_
        << nl;
    os  << indent << "multiply connected edges       : "
        << setw(8) << edges().size()-multipleStart_
        << nl;
    os  << decrIndent;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    Foam::extendedEdgeMesh::sideVolumeType& vt
)
{
    label type;
    is  >> type;

    vt = static_cast<Foam::extendedEdgeMesh::sideVolumeType>(type);

    // Check state of Istream
    is.check("operator>>(Istream&, sideVolumeType&)");

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Foam::extendedEdgeMesh::sideVolumeType& vt
)
{
    os  << static_cast<label>(vt);

    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const extendedEdgeMesh& em)
{
    // fileFormats::extendedEdgeMeshFormat::write(os, em.points_, em.edges_);
    os  << "// points" << nl
        << em.points() << nl
        << "// edges" << nl
        << em.edges() << nl
        << "// concaveStart mixedStart nonFeatureStart" << nl
        << em.concaveStart_ << token::SPACE
        << em.mixedStart_ << token::SPACE
        << em.nonFeatureStart_ << nl
        << "// internalStart flatStart openStart multipleStart" << nl
        << em.internalStart_ << token::SPACE
        << em.flatStart_ << token::SPACE
        << em.openStart_ << token::SPACE
        << em.multipleStart_ << nl
        << "// normals" << nl
        << em.normals_ << nl
        << "// normal volume types" << nl
        << em.normalVolumeTypes_ << nl
        << "// normalDirections" << nl
        << em.normalDirections_ << nl
        << "// edgeNormals" << nl
        << em.edgeNormals_ << nl
        << "// featurePointNormals" << nl
        << em.featurePointNormals_ << nl
        << "// featurePointEdges" << nl
        << em.featurePointEdges_ << nl
        << "// regionEdges" << nl
        << em.regionEdges_
        << endl;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const extendedEdgeMesh&)");

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, extendedEdgeMesh& em)
{
    // fileFormats::extendedEdgeMeshFormat::read(is, em.points_, em.edges_);
    is  >> static_cast<edgeMesh&>(em)
        >> em.concaveStart_
        >> em.mixedStart_
        >> em.nonFeatureStart_
        >> em.internalStart_
        >> em.flatStart_
        >> em.openStart_
        >> em.multipleStart_
        >> em.normals_
        >> em.normalVolumeTypes_
        >> em.normalDirections_
        >> em.edgeNormals_
        >> em.featurePointNormals_
        >> em.featurePointEdges_
        >> em.regionEdges_;

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, extendedEdgeMesh&)");

    return is;
}


// ************************************************************************* //
