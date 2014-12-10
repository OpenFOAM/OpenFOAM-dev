/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "dynamicIndexedOctree.H"
#include "linePointRef.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::dynamicIndexedOctree<Type>::perturbTol_ = 10*SMALL;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Does bb intersect a sphere around sample? Or is any corner point of bb
// closer than nearestDistSqr to sample.
template<class Type>
bool Foam::dynamicIndexedOctree<Type>::overlaps
(
    const point& p0,
    const point& p1,
    const scalar nearestDistSqr,
    const point& sample
)
{
    // Find out where sample is in relation to bb.
    // Find nearest point on bb.
    scalar distSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        scalar d0 = p0[dir] - sample[dir];
        scalar d1 = p1[dir] - sample[dir];

        if ((d0 > 0) != (d1 > 0))
        {
            // sample inside both extrema. This component does not add any
            // distance.
        }
        else if (mag(d0) < mag(d1))
        {
            distSqr += d0*d0;
        }
        else
        {
            distSqr += d1*d1;
        }

        if (distSqr > nearestDistSqr)
        {
            return false;
        }
    }

    return true;
}


// Does bb intersect a sphere around sample? Or is any corner point of bb
// closer than nearestDistSqr to sample.
template<class Type>
bool Foam::dynamicIndexedOctree<Type>::overlaps
(
    const treeBoundBox& parentBb,
    const direction octant,
    const scalar nearestDistSqr,
    const point& sample
)
{
    //- Accelerated version of
    //     treeBoundBox subBb(parentBb.subBbox(mid, octant))
    //     overlaps
    //     (
    //          subBb.min(),
    //          subBb.max(),
    //          nearestDistSqr,
    //          sample
    //     )

    const point& min = parentBb.min();
    const point& max = parentBb.max();

    point other;

    if (octant & treeBoundBox::RIGHTHALF)
    {
        other.x() = max.x();
    }
    else
    {
        other.x() = min.x();
    }

    if (octant & treeBoundBox::TOPHALF)
    {
        other.y() = max.y();
    }
    else
    {
        other.y() = min.y();
    }

    if (octant & treeBoundBox::FRONTHALF)
    {
        other.z() = max.z();
    }
    else
    {
        other.z() = min.z();
    }

    const point mid(0.5*(min+max));

    return overlaps(mid, other, nearestDistSqr, sample);
}


//
// Construction helper routines
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

// Split list of indices into 8 bins
template<class Type>
void Foam::dynamicIndexedOctree<Type>::divide
(
    const autoPtr<DynamicList<label> >& indices,
    const treeBoundBox& bb,
    contentListList& result
) const
{
    for (direction octant = 0; octant < 8; octant++)
    {
        result.append
        (
            autoPtr<DynamicList<label> >
            (
                new DynamicList<label>(indices().size()/8)
            )
        );
    }

    // Precalculate bounding boxes.
    FixedList<treeBoundBox, 8> subBbs;
    for (direction octant = 0; octant < 8; octant++)
    {
        subBbs[octant] = bb.subBbox(octant);
    }

    forAll(indices(), i)
    {
        label shapeI = indices()[i];

        for (direction octant = 0; octant < 8; octant++)
        {
            if (shapes_.overlaps(shapeI, subBbs[octant]))
            {
                result[octant]().append(shapeI);
            }
        }
    }
}


// Subdivide the (content) node.
template<class Type>
typename Foam::dynamicIndexedOctree<Type>::node
Foam::dynamicIndexedOctree<Type>::divide
(
    const treeBoundBox& bb,
    const label contentI,
    const label parentNodeIndex,
    const label octantToBeDivided
)
{
    const autoPtr<DynamicList<label> >& indices = contents_[contentI];

    node nod;

    if
    (
        bb.min()[0] >= bb.max()[0]
     || bb.min()[1] >= bb.max()[1]
     || bb.min()[2] >= bb.max()[2]
    )
    {
        FatalErrorIn("dynamicIndexedOctree<Type>::divide(..)")
            << "Badly formed bounding box:" << bb
            << abort(FatalError);
    }

    nod.bb_ = bb;
    nod.parent_ = -1;

    contentListList dividedIndices(8);
    divide(indices, bb, dividedIndices);

    // Have now divided the indices into 8 (possibly empty) subsets.
    // Replace current contentI with the first (non-empty) subset.
    // Append the rest.
    bool replaced = false;

    for (direction octant = 0; octant < dividedIndices.size(); octant++)
    {
        autoPtr<DynamicList<label> >& subIndices = dividedIndices[octant];

        if (subIndices().size())
        {
            if (!replaced)
            {
                contents_[contentI]().transfer(subIndices());
                nod.subNodes_[octant] = contentPlusOctant(contentI, octant);

                replaced = true;
            }
            else
            {
                // Store at end of contents.
                // note dummy append + transfer trick
                label sz = contents_.size();

                contents_.append
                (
                    autoPtr<DynamicList<label> >
                    (
                        new DynamicList<label>()
                    )
                );

                contents_[sz]().transfer(subIndices());

                nod.subNodes_[octant] = contentPlusOctant(sz, octant);
            }
        }
        else
        {
            // Mark node as empty
            nod.subNodes_[octant] = emptyPlusOctant(octant);
        }
    }

    // Don't update the parent node if it is the first node.
    if (parentNodeIndex != -1)
    {
        nod.parent_ = parentNodeIndex;

        label sz = nodes_.size();

        nodes_.append(nod);

        nodes_[parentNodeIndex].subNodes_[octantToBeDivided]
            = nodePlusOctant(sz, octantToBeDivided);
    }

    return nod;
}


template<class Type>
void Foam::dynamicIndexedOctree<Type>::recursiveSubDivision
(
    const treeBoundBox& subBb,
    const label contentI,
    const label parentIndex,
    const label octant,
    label& nLevels
)
{
    if
    (
        contents_[contentI]().size() > minSize_
     && nLevels < maxLevels_
    )
    {
        // Create node for content
        node nod = divide(subBb, contentI, parentIndex, octant);

        // Increment the number of levels in the tree
        nLevels++;

        // Recursively divide the contents until maxLevels_ is
        // reached or the content sizes are less than minSize_
        for (direction subOct = 0; subOct < 8; subOct++)
        {
            const labelBits& subNodeLabel = nod.subNodes_[subOct];

            if (isContent(subNodeLabel))
            {
                const treeBoundBox subBb = nod.bb_.subBbox(subOct);

                const label subContentI = getContent(subNodeLabel);

                const label parentNodeIndex = nodes_.size() - 1;

                recursiveSubDivision
                (
                    subBb,
                    subContentI,
                    parentNodeIndex,
                    subOct,
                    nLevels
                );
            }
        }
    }
}


// Pre-calculates wherever possible the volume status per node/subnode.
// Recurses to determine status of lowest level boxes. Level above is
// combination of octants below.
template<class Type>
Foam::volumeType Foam::dynamicIndexedOctree<Type>::calcVolumeType
(
    const label nodeI
) const
{
    const node& nod = nodes_[nodeI];

    volumeType myType = volumeType::UNKNOWN;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        volumeType subType;

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            // tree node. Recurse.
            subType = calcVolumeType(getNode(index));
        }
        else if (isContent(index))
        {
            // Contents. Depending on position in box might be on either
            // side.
            subType = volumeType::MIXED;
        }
        else
        {
            // No data in this octant. Set type for octant acc. to the mid
            // of its bounding box.
            const treeBoundBox subBb = nod.bb_.subBbox(octant);

            subType = volumeType
            (
                shapes_.getVolumeType(*this, subBb.midpoint())
            );
        }

        // Store octant type
        nodeTypes_.set((nodeI<<3)+octant, subType);

        // Combine sub node types into type for treeNode. Result is 'mixed' if
        // types differ among subnodes.
        if (myType == volumeType::UNKNOWN)
        {
            myType = subType;
        }
        else if (subType != myType)
        {
            myType = volumeType::MIXED;
        }
    }
    return myType;
}


template<class Type>
Foam::volumeType Foam::dynamicIndexedOctree<Type>::getVolumeType
(
    const label nodeI,
    const point& sample
) const
{
    const node& nod = nodes_[nodeI];

    direction octant = nod.bb_.subOctant(sample);

    volumeType octantType = volumeType::type(nodeTypes_.get((nodeI<<3)+octant));

    if (octantType == volumeType::INSIDE)
    {
        return octantType;
    }
    else if (octantType == volumeType::OUTSIDE)
    {
        return octantType;
    }
    else if (octantType == volumeType::UNKNOWN)
    {
        // Can happen for e.g. non-manifold surfaces.
        return octantType;
    }
    else if (octantType == volumeType::MIXED)
    {
        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            // Recurse
            volumeType subType = getVolumeType(getNode(index), sample);

            return subType;
        }
        else if (isContent(index))
        {
            // Content. Defer to shapes.
            return volumeType(shapes_.getVolumeType(*this, sample));
        }
        else
        {
            // Empty node. Cannot have 'mixed' as its type since not divided
            // up and has no items inside it.
            FatalErrorIn
            (
                "dynamicIndexedOctree<Type>::getVolumeType"
                "(const label, const point&)"
            )   << "Sample:" << sample << " node:" << nodeI
                << " with bb:" << nodes_[nodeI].bb_ << nl
                << "Empty subnode has invalid volume type MIXED."
                << abort(FatalError);

            return volumeType::UNKNOWN;
        }
    }
    else
    {
        FatalErrorIn
        (
            "dynamicIndexedOctree<Type>::getVolumeType"
            "(const label, const point&)"
        )   << "Sample:" << sample << " at node:" << nodeI
            << " octant:" << octant
            << " with bb:" << nod.bb_.subBbox(octant) << nl
            << "Node has invalid volume type " << octantType
            << abort(FatalError);

        return volumeType::UNKNOWN;
    }
}


template<class Type>
Foam::volumeType Foam::dynamicIndexedOctree<Type>::getSide
(
    const vector& outsideNormal,
    const vector& vec
)
{
    if ((outsideNormal&vec) >= 0)
    {
        return volumeType::OUTSIDE;
    }
    else
    {
        return volumeType::INSIDE;
    }
}


//
// Query routines
// ~~~~~~~~~~~~~~
//

// Find nearest point starting from nodeI
template<class Type>
void Foam::dynamicIndexedOctree<Type>::findNearest
(
    const label nodeI,
    const point& sample,

    scalar& nearestDistSqr,
    label& nearestShapeI,
    point& nearestPoint
) const
{
    const node& nod = nodes_[nodeI];

    // Determine order to walk through octants
    FixedList<direction, 8> octantOrder;
    nod.bb_.searchOrder(sample, octantOrder);

    // Go into all suboctants (one containing sample first) and update nearest.
    for (direction i = 0; i < 8; i++)
    {
        direction octant = octantOrder[i];

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            label subNodeI = getNode(index);

            const treeBoundBox& subBb = nodes_[subNodeI].bb_;

            if (overlaps(subBb.min(), subBb.max(), nearestDistSqr, sample))
            {
                findNearest
                (
                    subNodeI,
                    sample,

                    nearestDistSqr,
                    nearestShapeI,
                    nearestPoint
                );
            }
        }
        else if (isContent(index))
        {
            if
            (
                overlaps
                (
                    nod.bb_,
                    octant,
                    nearestDistSqr,
                    sample
                )
            )
            {
                shapes_.findNearest
                (
                    contents_[getContent(index)],
                    sample,

                    nearestDistSqr,
                    nearestShapeI,
                    nearestPoint
                );
            }
        }
    }
}


// Find nearest point to line.
template<class Type>
void Foam::dynamicIndexedOctree<Type>::findNearest
(
    const label nodeI,
    const linePointRef& ln,

    treeBoundBox& tightest,
    label& nearestShapeI,
    point& linePoint,
    point& nearestPoint
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& nodeBb = nod.bb_;

    // Determine order to walk through octants
    FixedList<direction, 8> octantOrder;
    nod.bb_.searchOrder(ln.centre(), octantOrder);

    // Go into all suboctants (one containing sample first) and update nearest.
    for (direction i = 0; i < 8; i++)
    {
        direction octant = octantOrder[i];

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            const treeBoundBox& subBb = nodes_[getNode(index)].bb_;

            if (subBb.overlaps(tightest))
            {
                findNearest
                (
                    getNode(index),
                    ln,

                    tightest,
                    nearestShapeI,
                    linePoint,
                    nearestPoint
                );
            }
        }
        else if (isContent(index))
        {
            const treeBoundBox subBb(nodeBb.subBbox(octant));

            if (subBb.overlaps(tightest))
            {
                shapes_.findNearest
                (
                    contents_[getContent(index)],
                    ln,

                    tightest,
                    nearestShapeI,
                    linePoint,
                    nearestPoint
                );
            }
        }
    }
}


template<class Type>
Foam::treeBoundBox Foam::dynamicIndexedOctree<Type>::subBbox
(
    const label parentNodeI,
    const direction octant
) const
{
    // Get type of node at octant
    const node& nod = nodes_[parentNodeI];
    labelBits index = nod.subNodes_[octant];

    if (isNode(index))
    {
        // Use stored bb
        return nodes_[getNode(index)].bb_;
    }
    else
    {
        // Calculate subBb
        return nod.bb_.subBbox(octant);
    }
}


// Takes a bb and a point on/close to the edge of the bb and pushes the point
// inside by a small fraction.
template<class Type>
Foam::point Foam::dynamicIndexedOctree<Type>::pushPoint
(
    const treeBoundBox& bb,
    const point& pt,
    const bool pushInside
)
{
    // Get local length scale.
    const vector perturbVec = perturbTol_*bb.span();

    point perturbedPt(pt);

    // Modify all components which are close to any face of the bb to be
    // well inside/outside them.

    if (pushInside)
    {
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (mag(pt[dir]-bb.min()[dir]) < mag(perturbVec[dir]))
            {
                // Close to 'left' side. Push well beyond left side.
                scalar perturbDist = perturbVec[dir] + ROOTVSMALL;
                perturbedPt[dir] = bb.min()[dir] + perturbDist;
            }
            else if (mag(pt[dir]-bb.max()[dir]) < mag(perturbVec[dir]))
            {
                // Close to 'right' side. Push well beyond right side.
                scalar perturbDist = perturbVec[dir] + ROOTVSMALL;
                perturbedPt[dir] = bb.max()[dir] - perturbDist;
            }
        }
    }
    else
    {
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (mag(pt[dir]-bb.min()[dir]) < mag(perturbVec[dir]))
            {
                scalar perturbDist = perturbVec[dir] + ROOTVSMALL;
                perturbedPt[dir] = bb.min()[dir] - perturbDist;
            }
            else if (mag(pt[dir]-bb.max()[dir]) < mag(perturbVec[dir]))
            {
                scalar perturbDist = perturbVec[dir] + ROOTVSMALL;
                perturbedPt[dir] = bb.max()[dir] + perturbDist;
            }
        }
    }

    if (debug)
    {
        if (pushInside != bb.contains(perturbedPt))
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::pushPoint(..)")
                << "pushed point:" << pt
                << " to:" << perturbedPt
                << " wanted side:" << pushInside
                << " obtained side:" << bb.contains(perturbedPt)
                << " of bb:" << bb
                << abort(FatalError);
        }
    }

    return perturbedPt;
}


// Takes a bb and a point on the edge of the bb and pushes the point
// outside by a small fraction.
template<class Type>
Foam::point Foam::dynamicIndexedOctree<Type>::pushPoint
(
    const treeBoundBox& bb,
    const direction faceID,
    const point& pt,
    const bool pushInside
)
{
    // Get local length scale.
    const vector perturbVec = perturbTol_*bb.span();

    point perturbedPt(pt);

    // Modify all components which are close to any face of the bb to be
    // well outside them.

    if (faceID == 0)
    {
        FatalErrorIn("dynamicIndexedOctree<Type>::pushPoint(..)")
            << abort(FatalError);
    }

    if (faceID & treeBoundBox::LEFTBIT)
    {
        if (pushInside)
        {
            perturbedPt[0] = bb.min()[0] + (perturbVec[0] + ROOTVSMALL);
        }
        else
        {
            perturbedPt[0] = bb.min()[0] - (perturbVec[0] + ROOTVSMALL);
        }
    }
    else if (faceID & treeBoundBox::RIGHTBIT)
    {
        if (pushInside)
        {
            perturbedPt[0] = bb.max()[0] - (perturbVec[0] + ROOTVSMALL);
        }
        else
        {
            perturbedPt[0] = bb.max()[0] + (perturbVec[0] + ROOTVSMALL);
        }
    }

    if (faceID & treeBoundBox::BOTTOMBIT)
    {
        if (pushInside)
        {
            perturbedPt[1] = bb.min()[1] + (perturbVec[1] + ROOTVSMALL);
        }
        else
        {
            perturbedPt[1] = bb.min()[1] - (perturbVec[1] + ROOTVSMALL);
        }
    }
    else if (faceID & treeBoundBox::TOPBIT)
    {
        if (pushInside)
        {
            perturbedPt[1] = bb.max()[1] - (perturbVec[1] + ROOTVSMALL);
        }
        else
        {
            perturbedPt[1] = bb.max()[1] + (perturbVec[1] + ROOTVSMALL);
        }
    }

    if (faceID & treeBoundBox::BACKBIT)
    {
        if (pushInside)
        {
            perturbedPt[2] = bb.min()[2] + (perturbVec[2] + ROOTVSMALL);
        }
        else
        {
            perturbedPt[2] = bb.min()[2] - (perturbVec[2] + ROOTVSMALL);
        }
    }
    else if (faceID & treeBoundBox::FRONTBIT)
    {
        if (pushInside)
        {
            perturbedPt[2] = bb.max()[2] - (perturbVec[2] + ROOTVSMALL);
        }
        else
        {
            perturbedPt[2] = bb.max()[2] + (perturbVec[2] + ROOTVSMALL);
        }
    }

    if (debug)
    {
        if (pushInside != bb.contains(perturbedPt))
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::pushPoint(..)")
                << "pushed point:" << pt << " on face:" << faceString(faceID)
                << " to:" << perturbedPt
                << " wanted side:" << pushInside
                << " obtained side:" << bb.contains(perturbedPt)
                << " of bb:" << bb
                << abort(FatalError);
        }
    }

    return perturbedPt;
}


// Guarantees that if pt is on a face it gets perturbed so it is away
// from the face edges.
// If pt is not on a face does nothing.
template<class Type>
Foam::point Foam::dynamicIndexedOctree<Type>::pushPointIntoFace
(
    const treeBoundBox& bb,
    const vector& dir,          // end-start
    const point& pt
)
{
    if (debug)
    {
        if (bb.posBits(pt) != 0)
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::pushPointIntoFace(..)")
                << " bb:" << bb << endl
                << "does not contain point " << pt << abort(FatalError);
        }
    }


    // Handle two cases:
    // - point exactly on multiple faces. Push away from all but one.
    // - point on a single face. Push away from edges of face.

    direction ptFaceID = bb.faceBits(pt);

    direction nFaces = 0;
    FixedList<direction, 3> faceIndices;

    if (ptFaceID & treeBoundBox::LEFTBIT)
    {
        faceIndices[nFaces++] = treeBoundBox::LEFT;
    }
    else if (ptFaceID & treeBoundBox::RIGHTBIT)
    {
        faceIndices[nFaces++] = treeBoundBox::RIGHT;
    }

    if (ptFaceID & treeBoundBox::BOTTOMBIT)
    {
        faceIndices[nFaces++] = treeBoundBox::BOTTOM;
    }
    else if (ptFaceID & treeBoundBox::TOPBIT)
    {
        faceIndices[nFaces++] = treeBoundBox::TOP;
    }

    if (ptFaceID & treeBoundBox::BACKBIT)
    {
        faceIndices[nFaces++] = treeBoundBox::BACK;
    }
    else if (ptFaceID & treeBoundBox::FRONTBIT)
    {
        faceIndices[nFaces++] = treeBoundBox::FRONT;
    }


    // Determine the face we want to keep the point on

    direction keepFaceID;

    if (nFaces == 0)
    {
        // Return original point
        return pt;
    }
    else if (nFaces == 1)
    {
        // Point is on a single face
        keepFaceID = faceIndices[0];
    }
    else
    {
        // Determine best face out of faceIndices[0 .. nFaces-1].
        // The best face is the one most perpendicular to the ray direction.

        keepFaceID = faceIndices[0];
        scalar maxInproduct = mag(treeBoundBox::faceNormals[keepFaceID] & dir);

        for (direction i = 1; i < nFaces; i++)
        {
            direction face = faceIndices[i];
            scalar s = mag(treeBoundBox::faceNormals[face] & dir);
            if (s > maxInproduct)
            {
                maxInproduct = s;
                keepFaceID = face;
            }
        }
    }


    // 1. Push point into bb, away from all corners

    point facePoint(pushPoint(bb, pt, true));
    direction faceID = 0;

    // 2. Snap it back onto the preferred face

    if (keepFaceID == treeBoundBox::LEFT)
    {
        facePoint.x() = bb.min().x();
        faceID = treeBoundBox::LEFTBIT;
    }
    else if (keepFaceID == treeBoundBox::RIGHT)
    {
        facePoint.x() = bb.max().x();
        faceID = treeBoundBox::RIGHTBIT;
    }
    else if (keepFaceID == treeBoundBox::BOTTOM)
    {
        facePoint.y() = bb.min().y();
        faceID = treeBoundBox::BOTTOMBIT;
    }
    else if (keepFaceID == treeBoundBox::TOP)
    {
        facePoint.y() = bb.max().y();
        faceID = treeBoundBox::TOPBIT;
    }
    else if (keepFaceID == treeBoundBox::BACK)
    {
        facePoint.z() = bb.min().z();
        faceID = treeBoundBox::BACKBIT;
    }
    else if (keepFaceID == treeBoundBox::FRONT)
    {
        facePoint.z() = bb.max().z();
        faceID = treeBoundBox::FRONTBIT;
    }


    if (debug)
    {
        if (faceID != bb.faceBits(facePoint))
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::pushPointIntoFace(..)")
                << "Pushed point from " << pt
                << " on face:" << ptFaceID << " of bb:" << bb << endl
                << "onto " << facePoint
                << " on face:" << faceID
                << " which is not consistent with geometric face "
                << bb.faceBits(facePoint)
                << abort(FatalError);
        }
        if (bb.posBits(facePoint) != 0)
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::pushPointIntoFace(..)")
                << " bb:" << bb << endl
                << "does not contain perturbed point "
                << facePoint << abort(FatalError);
        }
    }


    return facePoint;
}


//// Takes a bb and a point on the outside of the bb. Checks if on multiple
// faces
//// and if so perturbs point so it is only on one face.
//template<class Type>
//void Foam::dynamicIndexedOctree<Type>::checkMultipleFaces
//(
//    const treeBoundBox& bb,
//    const vector& dir,          // end-start
//    pointIndexHit& faceHitInfo,
//    direction& faceID
//)
//{
//    // Do the quick elimination of no or one face.
//    if
//    (
//        (faceID == 0)
//     || (faceID == treeBoundBox::LEFTBIT)
//     || (faceID == treeBoundBox::RIGHTBIT)
//     || (faceID == treeBoundBox::BOTTOMBIT)
//     || (faceID == treeBoundBox::TOPBIT)
//     || (faceID == treeBoundBox::BACKBIT)
//     || (faceID == treeBoundBox::FRONTBIT)
//    )
//    {
//        return;
//    }
//
//
//    // Check the direction of vector w.r.t. faces being intersected.
//    FixedList<scalar, 6> inproducts(-GREAT);
//
//    direction nFaces = 0;
//
//    if (faceID & treeBoundBox::LEFTBIT)
//    {
//        inproducts[treeBoundBox::LEFT] = mag
//        (
//            treeBoundBox::faceNormals[treeBoundBox::LEFT]
//          & dir
//        );
//        nFaces++;
//    }
//    if (faceID & treeBoundBox::RIGHTBIT)
//    {
//        inproducts[treeBoundBox::RIGHT] = mag
//        (
//            treeBoundBox::faceNormals[treeBoundBox::RIGHT]
//          & dir
//        );
//        nFaces++;
//    }
//
//    if (faceID & treeBoundBox::BOTTOMBIT)
//    {
//        inproducts[treeBoundBox::BOTTOM] = mag
//        (
//            treeBoundBox::faceNormals[treeBoundBox::BOTTOM]
//          & dir
//        );
//        nFaces++;
//    }
//    if (faceID & treeBoundBox::TOPBIT)
//    {
//        inproducts[treeBoundBox::TOP] = mag
//        (
//            treeBoundBox::faceNormals[treeBoundBox::TOP]
//          & dir
//        );
//        nFaces++;
//    }
//
//    if (faceID & treeBoundBox::BACKBIT)
//    {
//        inproducts[treeBoundBox::BACK] = mag
//        (
//            treeBoundBox::faceNormals[treeBoundBox::BACK]
//          & dir
//        );
//        nFaces++;
//    }
//    if (faceID & treeBoundBox::FRONTBIT)
//    {
//        inproducts[treeBoundBox::FRONT] = mag
//        (
//            treeBoundBox::faceNormals[treeBoundBox::FRONT]
//          & dir
//        );
//        nFaces++;
//    }
//
//    if (nFaces == 0 || nFaces == 1 || nFaces > 3)
//    {
//        FatalErrorIn("dynamicIndexedOctree<Type>::checkMultipleFaces(..)")
//            << "Problem : nFaces:" << nFaces << abort(FatalError);
//    }
//
//    // Keep point on most perpendicular face; shift it away from the aligned
//    // ones.
//    // E.g. line hits top and left face:
//    //     a
//    // ----+----+
//    //     |    |
//    //     |    |
//    //     +----+
//    // Shift point down (away from top):
//    //
//    //    a+----+
//    // ----|    |
//    //     |    |
//    //     +----+
//
//    label maxIndex = -1;
//    scalar maxInproduct = -GREAT;
//
//    for (direction i = 0; i < 6; i++)
//    {
//        if (inproducts[i] > maxInproduct)
//        {
//            maxInproduct = inproducts[i];
//            maxIndex = i;
//        }
//    }
//
//    if (maxIndex == -1)
//    {
//        FatalErrorIn("dynamicIndexedOctree<Type>::checkMultipleFaces(..)")
//            << "Problem maxIndex:" << maxIndex << " inproducts:" << inproducts
//            << abort(FatalError);
//    }
//
//    const point oldPoint(faceHitInfo.rawPoint());
//    const direction oldFaceID = faceID;
//
//    // 1. Push point into bb, away from all corners
//
//    faceHitInfo.rawPoint() = pushPoint(bb, oldFaceID, oldPoint, true);
//
//    // 2. Snap it back onto the preferred face
//
//    if (maxIndex == treeBoundBox::LEFT)
//    {
//        faceHitInfo.rawPoint().x() = bb.min().x();
//        faceID = treeBoundBox::LEFTBIT;
//    }
//    else if (maxIndex == treeBoundBox::RIGHT)
//    {
//        faceHitInfo.rawPoint().x() = bb.max().x();
//        faceID = treeBoundBox::RIGHTBIT;
//    }
//    else if (maxIndex == treeBoundBox::BOTTOM)
//    {
//        faceHitInfo.rawPoint().y() = bb.min().y();
//        faceID = treeBoundBox::BOTTOMBIT;
//    }
//    else if (maxIndex == treeBoundBox::TOP)
//    {
//        faceHitInfo.rawPoint().y() = bb.max().y();
//        faceID = treeBoundBox::TOPBIT;
//    }
//    else if (maxIndex == treeBoundBox::BACK)
//    {
//        faceHitInfo.rawPoint().z() = bb.min().z();
//        faceID = treeBoundBox::BACKBIT;
//    }
//    else if (maxIndex == treeBoundBox::FRONT)
//    {
//        faceHitInfo.rawPoint().z() = bb.max().z();
//        faceID = treeBoundBox::FRONTBIT;
//    }
//
//    Pout<< "From ray:" << dir
//        << " from point:" << oldPoint
//        << " on faces:" << faceString(oldFaceID)
//        << " of bb:" << bb
//        << " with inprods:" << inproducts
//        << " maxIndex:" << maxIndex << endl
//        << "perturbed to point:" << faceHitInfo.rawPoint()
//        << " on face:" << faceString(faceID)
//        << endl;
//
//
//    if (debug)
//    {
//        if (faceID != bb.faceBits(faceHitInfo.rawPoint()))
//        {
//            FatalErrorIn("dynamicIndexedOctree<Type>::checkMultipleFaces(..)")
//                << "Pushed point from " << oldPoint
//                << " on face:" << oldFaceID << " of bb:" << bb << endl
//                << "onto " << faceHitInfo.rawPoint()
//                << " on face:" << faceID
//                << " which is not consistent with geometric face "
//                <<  bb.faceBits(faceHitInfo.rawPoint())
//                << abort(FatalError);
//        }
//    }
//}


// Get parent node and octant. Return false if top of tree reached.
template<class Type>
bool Foam::dynamicIndexedOctree<Type>::walkToParent
(
    const label nodeI,
    const direction octant,

    label& parentNodeI,
    label& parentOctant
) const
{
    parentNodeI = nodes_[nodeI].parent_;

    if (parentNodeI == -1)
    {
        // Reached edge of tree
        return false;
    }

    const node& parentNode = nodes_[parentNodeI];

    // Find octant nodeI is in.
    parentOctant = 255;

    for (direction i = 0; i < parentNode.subNodes_.size(); i++)
    {
        labelBits index = parentNode.subNodes_[i];

        if (isNode(index) && getNode(index) == nodeI)
        {
            parentOctant = i;
            break;
        }
    }

    if (parentOctant == 255)
    {
        FatalErrorIn("walkToParent(..)")
            << "Problem: no parent found for octant:" << octant
            << " node:" << nodeI
            << abort(FatalError);
    }

    return true;
}


// Walk tree to neighbouring node. Gets current position as
// node and octant in this node and walks in the direction given by
// the facePointBits (combination of treeBoundBox::LEFTBIT, TOPBIT etc.)
// Returns false if edge of tree hit.
template<class Type>
bool Foam::dynamicIndexedOctree<Type>::walkToNeighbour
(
    const point& facePoint,
    const direction faceID,  // face(s) that facePoint is on
    label& nodeI,
    direction& octant
) const
{
    label oldNodeI = nodeI;
    direction oldOctant = octant;

    // Find out how to test for octant. Say if we want to go left we need
    // to traverse up the tree until we hit a node where our octant is
    // on the right.

    // Coordinate direction to test
    const direction X = treeBoundBox::RIGHTHALF;
    const direction Y = treeBoundBox::TOPHALF;
    const direction Z = treeBoundBox::FRONTHALF;

    direction octantMask = 0;
    direction wantedValue = 0;

    if ((faceID & treeBoundBox::LEFTBIT) != 0)
    {
        // We want to go left so check if in right octant (i.e. x-bit is set)
        octantMask |= X;
        wantedValue |= X;
    }
    else if ((faceID & treeBoundBox::RIGHTBIT) != 0)
    {
        octantMask |= X;  // wantedValue already 0
    }

    if ((faceID & treeBoundBox::BOTTOMBIT) != 0)
    {
        // Want to go down so check for y-bit set.
        octantMask |= Y;
        wantedValue |= Y;
    }
    else if ((faceID & treeBoundBox::TOPBIT) != 0)
    {
        // Want to go up so check for y-bit not set.
        octantMask |= Y;
    }

    if ((faceID & treeBoundBox::BACKBIT) != 0)
    {
        octantMask |= Z;
        wantedValue |= Z;
    }
    else if ((faceID & treeBoundBox::FRONTBIT) != 0)
    {
        octantMask |= Z;
    }

    // So now we have the coordinate directions in the octant we need to check
    // and the resulting value.

    /*
    // +---+---+
    // |   |   |
    // |   |   |
    // |   |   |
    // +---+-+-+
    // |   | | |
    // |  a+-+-+
    // |   |\| |
    // +---+-+-+
    //        \
    //
    // e.g. ray is at (a) in octant 0(or 4) with faceIDs : LEFTBIT+TOPBIT.
    // If we would be in octant 1(or 5) we could go to the correct octant
    // in the same node by just flipping the x and y bits (exoring).
    // But if we are not in octant 1/5 we have to go up until we are.
    // In general for leftbit+topbit:
    // - we have to check for x and y : octantMask  = 011
    // - the checked bits have to be  : wantedValue = ?01
    */

    //Pout<< "For point " << facePoint << endl;

    // Go up until we have chance to cross to the wanted direction
    while (wantedValue != (octant & octantMask))
    {
        // Go up to the parent.

        // Remove the directions that are not on the boundary of the parent.
        // See diagram above
        if (wantedValue & X)    // && octantMask&X
        {
            // Looking for right octant.
            if (octant & X)
            {
                // My octant is one of the left ones so punch out the x check
                octantMask &= ~X;
                wantedValue &= ~X;
            }
        }
        else
        {
            if (!(octant & X))
            {
                // My octant is right but I want to go left.
                octantMask &= ~X;
                wantedValue &= ~X;
            }
        }

        if (wantedValue & Y)
        {
            if (octant & Y)
            {
                octantMask &= ~Y;
                wantedValue &= ~Y;
            }
        }
        else
        {
            if (!(octant & Y))
            {
                octantMask &= ~Y;
                wantedValue &= ~Y;
            }
        }

        if (wantedValue & Z)
        {
            if (octant & Z)
            {
                octantMask &= ~Z;
                wantedValue &= ~Z;
            }
        }
        else
        {
            if (!(octant & Z))
            {
                octantMask &= ~Z;
                wantedValue &= ~Z;
            }
        }


        label parentNodeI;
        label parentOctant;
        walkToParent(nodeI, octant, parentNodeI, parentOctant);

        if (parentNodeI == -1)
        {
            // Reached edge of tree
            return false;
        }

        //Pout<< "    walked from node:" << nodeI << " octant:" << octant
        //    << " bb:" << nodes_[nodeI].bb_.subBbox(octant) << endl
        //    << "    to:" << parentNodeI << " octant:" << parentOctant
        //    << " bb:" << nodes_[parentNodeI].bb_.subBbox(parentOctant)
        //    << endl;
        //
        //Pout<< "    octantMask:" << octantMask
        //    << " wantedValue:" << wantedValue << endl;

        nodeI = parentNodeI;
        octant = parentOctant;
    }

    // So now we hit the correct parent node. Switch to the correct octant.
    // Is done by jumping to the 'other half' so if e.g. in x direction in
    // right half we now jump to the left half.
    octant ^= octantMask;

    //Pout<< "    to node:" << nodeI << " octant:" << octant
    //    << " subBb:" <<subBbox(nodeI, octant) << endl;


    if (debug)
    {
        const treeBoundBox subBb(subBbox(nodeI, octant));

        if (!subBb.contains(facePoint))
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::walkToNeighbour(..)")
                << "When searching for " << facePoint
                << " ended up in node:" << nodeI
                << " octant:" << octant
                << " with bb:" << subBb
                << abort(FatalError);
        }
    }


    // See if we need to travel down. Note that we already go into the
    // the first level ourselves (instead of having findNode decide)
    labelBits index = nodes_[nodeI].subNodes_[octant];

    if (isNode(index))
    {
        labelBits node = findNode(getNode(index), facePoint);

        nodeI = getNode(node);
        octant = getOctant(node);
    }


    if (debug)
    {
        const treeBoundBox subBb(subBbox(nodeI, octant));

        if (nodeI == oldNodeI && octant == oldOctant)
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::walkToNeighbour(..)")
                << "Did not go to neighbour when searching for " << facePoint
                << endl
                << "    starting from face:" << faceString(faceID)
                << " node:" << nodeI
                << " octant:" << octant
                << " bb:" << subBb
                << abort(FatalError);
        }

        if (!subBb.contains(facePoint))
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::walkToNeighbour(..)")
                << "When searching for " << facePoint
                << " ended up in node:" << nodeI
                << " octant:" << octant
                << " bb:" << subBb
                << abort(FatalError);
        }
    }


    return true;
}


template<class Type>
Foam::word Foam::dynamicIndexedOctree<Type>::faceString
(
    const direction faceID
)
{
    word desc;

    if (faceID == 0)
    {
        desc = "noFace";
    }
    if (faceID & treeBoundBox::LEFTBIT)
    {
        if (!desc.empty()) desc += "+";
        desc += "left";
    }
    if (faceID & treeBoundBox::RIGHTBIT)
    {
        if (!desc.empty()) desc += "+";
        desc += "right";
    }
    if (faceID & treeBoundBox::BOTTOMBIT)
    {
        if (!desc.empty()) desc += "+";
        desc += "bottom";
    }
    if (faceID & treeBoundBox::TOPBIT)
    {
        if (!desc.empty()) desc += "+";
        desc += "top";
    }
    if (faceID & treeBoundBox::BACKBIT)
    {
        if (!desc.empty()) desc += "+";
        desc += "back";
    }
    if (faceID & treeBoundBox::FRONTBIT)
    {
        if (!desc.empty()) desc += "+";
        desc += "front";
    }
    return desc;
}


// Traverse a node. If intersects a triangle return first intersection point:
//  hitInfo.index = index of shape
//  hitInfo.point = point on shape
// Else return a miss and the bounding box face hit:
//  hitInfo.point = coordinate of intersection of ray with bounding box
//  hitBits  = posbits of point on bounding box
template<class Type>
void Foam::dynamicIndexedOctree<Type>::traverseNode
(
    const bool findAny,
    const point& treeStart,
    const vector& treeVec,

    const point& start,
    const point& end,
    const label nodeI,
    const direction octant,

    pointIndexHit& hitInfo,
    direction& hitBits
) const
{
    if (debug)
    {
        const treeBoundBox octantBb(subBbox(nodeI, octant));

        if (octantBb.posBits(start) != 0)
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::traverseNode(..)")
                << "Node:" << nodeI << " octant:" << octant
                << " bb:" << octantBb << endl
                << "does not contain point " << start << abort(FatalError);
        }
    }


    const node& nod = nodes_[nodeI];

    labelBits index = nod.subNodes_[octant];

    if (isContent(index))
    {
        const labelList& indices = contents_[getContent(index)];

        if (indices.size())
        {
            if (findAny)
            {
                // Find any intersection

                forAll(indices, elemI)
                {
                    label shapeI = indices[elemI];

                    point pt;
                    bool hit = shapes_.intersects(shapeI, start, end, pt);

                    // Note that intersection of shape might actually be
                    // in a neighbouring box. For findAny this is not important.
                    if (hit)
                    {
                        // Hit so pt is nearer than nearestPoint.
                        // Update hit info
                        hitInfo.setHit();
                        hitInfo.setIndex(shapeI);
                        hitInfo.setPoint(pt);
                        return;
                    }
                }
            }
            else
            {
                // Find nearest intersection

                const treeBoundBox octantBb(subBbox(nodeI, octant));

                point nearestPoint(end);

                forAll(indices, elemI)
                {
                    label shapeI = indices[elemI];

                    point pt;
                    bool hit = shapes_.intersects
                    (
                        shapeI,
                        start,
                        nearestPoint,
                        pt
                    );

                    // Note that intersection of shape might actually be
                    // in a neighbouring box. Since we need to maintain strict
                    // (findAny=false) ordering skip such an intersection. It
                    // will be found when we are doing the next box.

                    if (hit && octantBb.contains(pt))
                    {
                        // Hit so pt is nearer than nearestPoint.
                        nearestPoint = pt;
                        // Update hit info
                        hitInfo.setHit();
                        hitInfo.setIndex(shapeI);
                        hitInfo.setPoint(pt);
                    }
                }

                if (hitInfo.hit())
                {
                    // Found intersected shape.
                    return;
                }
            }
        }
    }

    // Nothing intersected in this node
    // Traverse to end of node. Do by ray tracing back from end and finding
    // intersection with bounding box of node.
    // start point is now considered to be inside bounding box.

    const treeBoundBox octantBb(subBbox(nodeI, octant));

    point pt;
    bool intersected = octantBb.intersects
    (
        end,            //treeStart,
        (start-end),    //treeVec,

        end,
        start,

        pt,
        hitBits
    );


    if (intersected)
    {
        // Return miss. Set misspoint to face.
        hitInfo.setPoint(pt);
    }
    else
    {
        // Rare case: did not find intersection of ray with octantBb. Can
        // happen if end is on face/edge of octantBb. Do like
        // lagrangian and re-track with perturbed vector (slightly pushed out
        // of bounding box)

        point perturbedEnd(pushPoint(octantBb, end, false));

        traverseNode
        (
            findAny,
            treeStart,
            treeVec,
            start,
            perturbedEnd,
            nodeI,
            octant,

            hitInfo,
            hitBits
        );
    }
}


// Find first intersection
template<class Type>
Foam::pointIndexHit Foam::dynamicIndexedOctree<Type>::findLine
(
    const bool findAny,
    const point& treeStart,
    const point& treeEnd,
    const label startNodeI,
    const direction startOctant,
    const bool verbose
) const
{
    const vector treeVec(treeEnd - treeStart);

    // Current node as parent+octant
    label nodeI = startNodeI;
    direction octant = startOctant;

    if (verbose)
    {
        Pout<< "findLine : treeStart:" << treeStart
            << " treeEnd:" << treeEnd << endl
            << "node:" << nodeI
            << " octant:" << octant
            << " bb:" << subBbox(nodeI, octant) << endl;
    }

    // Current position. Initialize to miss
    pointIndexHit hitInfo(false, treeStart, -1);

    //while (true)
    label i = 0;
    for (; i < 100000; i++)
    {
        // Ray-trace to end of current node. Updates point (either on triangle
        // in case of hit or on node bounding box in case of miss)

        const treeBoundBox octantBb(subBbox(nodeI, octant));

        // Make sure point is away from any edges/corners
        point startPoint
        (
            pushPointIntoFace
            (
                octantBb,
                treeVec,
                hitInfo.rawPoint()
            )
        );

        if (verbose)
        {
            Pout<< "iter:" << i
                << " at current:" << hitInfo.rawPoint()
                << " (perturbed:" << startPoint << ")" << endl
                << "    node:" << nodeI
                << " octant:" << octant
                << " bb:" << subBbox(nodeI, octant) << endl;
        }




        // Faces of current bounding box current point is on
        direction hitFaceID = 0;

        traverseNode
        (
            findAny,
            treeStart,
            treeVec,

            startPoint,     // Note: pass in copy since hitInfo
                            // also used in return value.
            treeEnd,        // pass in overall end so is nicely outside bb
            nodeI,
            octant,

            hitInfo,
            hitFaceID
        );

        // Did we hit a triangle?
        if (hitInfo.hit())
        {
            break;
        }

        if (hitFaceID == 0 || hitInfo.rawPoint() == treeEnd)
        {
            // endpoint inside the tree. Return miss.
            break;
        }

        // Create a point on other side of face.
        point perturbedPoint
        (
            pushPoint
            (
                octantBb,
                hitFaceID,
                hitInfo.rawPoint(),
                false                   // push outside of octantBb
            )
        );

        if (verbose)
        {
            Pout<< "    iter:" << i
                << " hit face:" << faceString(hitFaceID)
                << " at:" << hitInfo.rawPoint() << nl
                << "    node:" << nodeI
                << " octant:" << octant
                << " bb:" << subBbox(nodeI, octant) << nl
                << "    walking to neighbour containing:" << perturbedPoint
                << endl;
        }


        // Nothing hit so we are on face of bounding box (given as node+octant+
        // position bits). Traverse to neighbouring node. Use slightly perturbed
        // point.

        bool ok = walkToNeighbour
        (
            perturbedPoint,
            hitFaceID,  // face(s) that hitInfo is on

            nodeI,
            octant
        );

        if (!ok)
        {
            // Hit the edge of the tree. Return miss.
            break;
        }

        if (verbose)
        {
            const treeBoundBox octantBb(subBbox(nodeI, octant));
            Pout<< "    walked for point:" << hitInfo.rawPoint() << endl
                << "    to neighbour node:" << nodeI
                << " octant:" << octant
                << " face:" << faceString(octantBb.faceBits(hitInfo.rawPoint()))
                << " of octantBb:" << octantBb << endl
                << endl;
        }
    }

    if (i == 100000)
    {
        // Probably in loop.
        if (!verbose)
        {
            // Redo intersection but now with verbose flag switched on.
            return findLine
            (
                findAny,
                treeStart,
                treeEnd,
                startNodeI,
                startOctant,
                true            //verbose
            );
        }
        if (debug)
        {
            FatalErrorIn("dynamicIndexedOctree<Type>::findLine(..)")
                << "Got stuck in loop raytracing from:" << treeStart
                << " to:" << treeEnd << endl
                << "inside top box:" << subBbox(startNodeI, startOctant)
                << abort(FatalError);
        }
        else
        {
            WarningIn("dynamicIndexedOctree<Type>::findLine(..)")
                << "Got stuck in loop raytracing from:" << treeStart
                << " to:" << treeEnd << endl
                << "inside top box:" << subBbox(startNodeI, startOctant)
                << endl;
        }
    }

    return hitInfo;
}


// Find first intersection
template<class Type>
Foam::pointIndexHit Foam::dynamicIndexedOctree<Type>::findLine
(
    const bool findAny,
    const point& start,
    const point& end
) const
{
    pointIndexHit hitInfo;

    if (nodes_.size())
    {
        const treeBoundBox& treeBb = nodes_[0].bb_;

        // No effort is made to deal with points which are on edge of tree
        // bounding box for now.

        direction startBit = treeBb.posBits(start);
        direction endBit = treeBb.posBits(end);

        if ((startBit & endBit) != 0)
        {
            // Both start and end outside domain and in same block.
            return pointIndexHit(false, vector::zero, -1);
        }


        // Trim segment to treeBb

        point trackStart(start);
        point trackEnd(end);

        if (startBit != 0)
        {
            // Track start to inside domain.
            if (!treeBb.intersects(start, end, trackStart))
            {
                return pointIndexHit(false, vector::zero, -1);
            }
        }

        if (endBit != 0)
        {
            // Track end to inside domain.
            if (!treeBb.intersects(end, trackStart, trackEnd))
            {
                return pointIndexHit(false, vector::zero, -1);
            }
        }


        // Find lowest level tree node that start is in.
        labelBits index = findNode(0, trackStart);

        label parentNodeI = getNode(index);
        direction octant = getOctant(index);

        hitInfo = findLine
        (
            findAny,
            trackStart,
            trackEnd,
            parentNodeI,
            octant
        );
    }

    return hitInfo;
}


template<class Type>
void Foam::dynamicIndexedOctree<Type>::findBox
(
    const label nodeI,
    const treeBoundBox& searchBox,
    labelHashSet& elements
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& nodeBb = nod.bb_;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            const treeBoundBox& subBb = nodes_[getNode(index)].bb_;

            if (subBb.overlaps(searchBox))
            {
                findBox(getNode(index), searchBox, elements);
            }
        }
        else if (isContent(index))
        {
            const treeBoundBox subBb(nodeBb.subBbox(octant));

            if (subBb.overlaps(searchBox))
            {
                const labelList& indices = contents_[getContent(index)];

                forAll(indices, i)
                {
                    label shapeI = indices[i];

                    if (shapes_.overlaps(shapeI, searchBox))
                    {
                        elements.insert(shapeI);
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::dynamicIndexedOctree<Type>::findSphere
(
    const label nodeI,
    const point& centre,
    const scalar radiusSqr,
    labelHashSet& elements
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& nodeBb = nod.bb_;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            const treeBoundBox& subBb = nodes_[getNode(index)].bb_;

            if (subBb.overlaps(centre, radiusSqr))
            {
                findSphere(getNode(index), centre, radiusSqr, elements);
            }
        }
        else if (isContent(index))
        {
            const treeBoundBox subBb(nodeBb.subBbox(octant));

            if (subBb.overlaps(centre, radiusSqr))
            {
                const labelList& indices = contents_[getContent(index)];

                forAll(indices, i)
                {
                    label shapeI = indices[i];

                    if (shapes_.overlaps(shapeI, centre, radiusSqr))
                    {
                        elements.insert(shapeI);
                    }
                }
            }
        }
    }
}


template<class Type>
template<class CompareOp>
void Foam::dynamicIndexedOctree<Type>::findNear
(
    const scalar nearDist,
    const bool okOrder,
    const dynamicIndexedOctree<Type>& tree1,
    const labelBits index1,
    const treeBoundBox& bb1,
    const dynamicIndexedOctree<Type>& tree2,
    const labelBits index2,
    const treeBoundBox& bb2,
    CompareOp& cop
)
{
    const vector nearSpan(nearDist, nearDist, nearDist);

    if (tree1.isNode(index1))
    {
        const node& nod1 = tree1.nodes()[tree1.getNode(index1)];
        const treeBoundBox searchBox
        (
            bb1.min()-nearSpan,
            bb1.max()+nearSpan
        );

        if (tree2.isNode(index2))
        {
            if (bb2.overlaps(searchBox))
            {
                const node& nod2 = tree2.nodes()[tree2.getNode(index2)];

                for (direction i2 = 0; i2 < nod2.subNodes_.size(); i2++)
                {
                    labelBits subIndex2 = nod2.subNodes_[i2];
                    const treeBoundBox subBb2
                    (
                        tree2.isNode(subIndex2)
                      ? tree2.nodes()[tree2.getNode(subIndex2)].bb_
                      : bb2.subBbox(i2)
                    );

                    findNear
                    (
                        nearDist,
                        !okOrder,
                        tree2,
                        subIndex2,
                        subBb2,
                        tree1,
                        index1,
                        bb1,
                        cop
                    );
                }
            }
        }
        else if (tree2.isContent(index2))
        {
            // index2 is leaf, index1 not yet.
            for (direction i1 = 0; i1 < nod1.subNodes_.size(); i1++)
            {
                labelBits subIndex1 = nod1.subNodes_[i1];
                const treeBoundBox subBb1
                (
                    tree1.isNode(subIndex1)
                  ? tree1.nodes()[tree1.getNode(subIndex1)].bb_
                  : bb1.subBbox(i1)
                );

                findNear
                (
                    nearDist,
                    !okOrder,
                    tree2,
                    index2,
                    bb2,
                    tree1,
                    subIndex1,
                    subBb1,
                    cop
                );
            }
        }
    }
    else if (tree1.isContent(index1))
    {
        const treeBoundBox searchBox
        (
            bb1.min()-nearSpan,
            bb1.max()+nearSpan
        );

        if (tree2.isNode(index2))
        {
            const node& nod2 =
                tree2.nodes()[tree2.getNode(index2)];

            if (bb2.overlaps(searchBox))
            {
                for (direction i2 = 0; i2 < nod2.subNodes_.size(); i2++)
                {
                    labelBits subIndex2 = nod2.subNodes_[i2];
                    const treeBoundBox subBb2
                    (
                        tree2.isNode(subIndex2)
                      ? tree2.nodes()[tree2.getNode(subIndex2)].bb_
                      : bb2.subBbox(i2)
                    );

                    findNear
                    (
                        nearDist,
                        !okOrder,
                        tree2,
                        subIndex2,
                        subBb2,
                        tree1,
                        index1,
                        bb1,
                        cop
                    );
                }
            }
        }
        else if (tree2.isContent(index2))
        {
            // Both are leaves. Check n^2.

            const labelList& indices1 =
                tree1.contents()[tree1.getContent(index1)];
            const labelList& indices2 =
                tree2.contents()[tree2.getContent(index2)];

            forAll(indices1, i)
            {
                label shape1 = indices1[i];

                forAll(indices2, j)
                {
                    label shape2 = indices2[j];

                    if ((&tree1 != &tree2) || (shape1 != shape2))
                    {
                        if (okOrder)
                        {
                            cop
                            (
                                nearDist,
                                tree1.shapes(),
                                shape1,
                                tree2.shapes(),
                                shape2
                            );
                        }
                        else
                        {
                            cop
                            (
                                nearDist,
                                tree2.shapes(),
                                shape2,
                                tree1.shapes(),
                                shape1
                            );
                        }
                    }
                }
            }
        }
    }
}


// Number of elements in node.
template<class Type>
Foam::label Foam::dynamicIndexedOctree<Type>::countElements
(
    const labelBits index
) const
{
    label nElems = 0;

    if (isNode(index))
    {
        // tree node.
        label nodeI = getNode(index);

        const node& nod = nodes_[nodeI];

        for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
        {
            nElems += countElements(nod.subNodes_[octant]);
        }
    }
    else if (isContent(index))
    {
        nElems += contents_[getContent(index)]().size();
    }
    else
    {
        // empty node
    }

    return nElems;
}


template<class Type>
void Foam::dynamicIndexedOctree<Type>::writeOBJ
(
    const label nodeI,
    const direction octant
) const
{
    OFstream str
    (
        "node" + Foam::name(nodeI) + "_octant" + Foam::name(octant) + ".obj"
    );

    labelBits index = nodes_[nodeI].subNodes_[octant];

    treeBoundBox subBb;

    if (isNode(index))
    {
        subBb = nodes_[getNode(index)].bb_;
    }
    else if (isContent(index) || isEmpty(index))
    {
        subBb = nodes_[nodeI].bb_.subBbox(octant);
    }

    Pout<< "dumpContentNode : writing node:" << nodeI << " octant:" << octant
        << " to " << str.name() << endl;

    // Dump bounding box
    pointField bbPoints(subBb.points());

    forAll(bbPoints, i)
    {
        const point& pt = bbPoints[i];

        str<< "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0] + 1 << ' ' << e[1] + 1 << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::dynamicIndexedOctree<Type>::dynamicIndexedOctree
(
    const Type& shapes,
    const treeBoundBox& bb,
    const label maxLevels,          // maximum number of levels
    const scalar maxLeafRatio,
    const scalar maxDuplicity
)
:
    shapes_(shapes),
    bb_(bb),
    maxLevels_(maxLevels),
    nLevelsMax_(0),
    maxLeafRatio_(maxLeafRatio),
    minSize_(label(maxLeafRatio)),
    maxDuplicity_(maxDuplicity),
    nodes_(label(shapes.size() / maxLeafRatio_)),
    contents_(label(shapes.size() / maxLeafRatio_)),
    nodeTypes_(0)
{
    if (shapes_.size() == 0)
    {
        return;
    }

    insert(0, shapes_.size());

    if (debug)
    {
        writeTreeInfo();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::scalar& Foam::dynamicIndexedOctree<Type>::perturbTol()
{
    return perturbTol_;
}


template<class Type>
Foam::pointIndexHit Foam::dynamicIndexedOctree<Type>::findNearest
(
    const point& sample,
    const scalar startDistSqr
) const
{
    scalar nearestDistSqr = startDistSqr;
    label nearestShapeI = -1;
    point nearestPoint = vector::zero;

    if (nodes_.size())
    {
        findNearest
        (
            0,
            sample,

            nearestDistSqr,
            nearestShapeI,
            nearestPoint
        );
    }

    return pointIndexHit(nearestShapeI != -1, nearestPoint, nearestShapeI);
}


template<class Type>
Foam::pointIndexHit Foam::dynamicIndexedOctree<Type>::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    label nearestShapeI = -1;
    point nearestPoint;

    if (nodes_.size())
    {
        findNearest
        (
            0,
            ln,

            tightest,
            nearestShapeI,
            linePoint,
            nearestPoint
        );
    }
    else
    {
        nearestPoint = vector::zero;
    }

    return pointIndexHit(nearestShapeI != -1, nearestPoint, nearestShapeI);
}


// Find nearest intersection
template<class Type>
Foam::pointIndexHit Foam::dynamicIndexedOctree<Type>::findLine
(
    const point& start,
    const point& end
) const
{
    return findLine(false, start, end);
}


// Find nearest intersection
template<class Type>
Foam::pointIndexHit Foam::dynamicIndexedOctree<Type>::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine(true, start, end);
}


template<class Type>
Foam::labelList Foam::dynamicIndexedOctree<Type>::findBox
(
    const treeBoundBox& searchBox
) const
{
    // Storage for labels of shapes inside bb. Size estimate.
    labelHashSet elements(shapes_.size() / 100);

    if (nodes_.size())
    {
        findBox(0, searchBox, elements);
    }

    return elements.toc();
}


template<class Type>
Foam::labelList Foam::dynamicIndexedOctree<Type>::findSphere
(
    const point& centre,
    const scalar radiusSqr
) const
{
    // Storage for labels of shapes inside bb. Size estimate.
    labelHashSet elements(shapes_.size() / 100);

    if (nodes_.size())
    {
        findSphere(0, centre, radiusSqr, elements);
    }

    return elements.toc();
}


// Find node (as parent+octant) containing point
template<class Type>
Foam::labelBits Foam::dynamicIndexedOctree<Type>::findNode
(
    const label nodeI,
    const point& sample
) const
{
    if (nodes_.empty())
    {
        // Empty tree. Return what?
        return nodePlusOctant(nodeI, 0);
    }

    const node& nod = nodes_[nodeI];

    if (debug)
    {
        if (!nod.bb_.contains(sample))
        {
            FatalErrorIn("findNode(..)")
                << "Cannot find " << sample << " in node " << nodeI
                << abort(FatalError);
        }
    }

    direction octant = nod.bb_.subOctant(sample);

    labelBits index = nod.subNodes_[octant];

    if (isNode(index))
    {
        // Recurse
        return findNode(getNode(index), sample);
    }
    else if (isContent(index))
    {
        // Content. Return treenode+octant
        return nodePlusOctant(nodeI, octant);
    }
    else
    {
        // Empty. Return treenode+octant
        return nodePlusOctant(nodeI, octant);
    }
}


template<class Type>
Foam::label Foam::dynamicIndexedOctree<Type>::findInside
(
    const point& sample
) const
{
    labelBits index = findNode(0, sample);

    const node& nod = nodes_[getNode(index)];

    labelBits contentIndex = nod.subNodes_[getOctant(index)];

    // Need to check for the presence of content, in-case the node is empty
    if (isContent(contentIndex))
    {
        labelList indices = contents_[getContent(contentIndex)];

        forAll(indices, elemI)
        {
            label shapeI = indices[elemI];

            if (shapes_.contains(shapeI, sample))
            {
                return shapeI;
            }
        }
    }

    return -1;
}


template<class Type>
const Foam::labelList& Foam::dynamicIndexedOctree<Type>::findIndices
(
    const point& sample
) const
{
    labelBits index = findNode(0, sample);

    const node& nod = nodes_[getNode(index)];

    labelBits contentIndex = nod.subNodes_[getOctant(index)];

    // Need to check for the presence of content, in-case the node is empty
    if (isContent(contentIndex))
    {
        return contents_[getContent(contentIndex)];
    }
    else
    {
        return emptyList<label>();
    }
}


// Determine type (inside/outside/mixed) per node.
template<class Type>
Foam::volumeType Foam::dynamicIndexedOctree<Type>::getVolumeType
(
    const point& sample
) const
{
    if (nodes_.empty())
    {
        return volumeType::UNKNOWN;
    }

    if (nodeTypes_.size() != 8*nodes_.size())
    {
        // Calculate type for every octant of node.

        nodeTypes_.setSize(8*nodes_.size());
        nodeTypes_ = volumeType::UNKNOWN;

        calcVolumeType(0);

        if (debug)
        {
            label nUNKNOWN = 0;
            label nMIXED = 0;
            label nINSIDE = 0;
            label nOUTSIDE = 0;

            forAll(nodeTypes_, i)
            {
                volumeType type = volumeType::type(nodeTypes_.get(i));

                if (type == volumeType::UNKNOWN)
                {
                    nUNKNOWN++;
                }
                else if (type == volumeType::MIXED)
                {
                    nMIXED++;
                }
                else if (type == volumeType::INSIDE)
                {
                    nINSIDE++;
                }
                else if (type == volumeType::OUTSIDE)
                {
                    nOUTSIDE++;
                }
                else
                {
                    FatalErrorIn("getVolumeType") << abort(FatalError);
                }
            }

            Pout<< "dynamicIndexedOctree<Type>::getVolumeType : "
                << " bb:" << bb()
                << " nodes_:" << nodes_.size()
                << " nodeTypes_:" << nodeTypes_.size()
                << " nUNKNOWN:" << nUNKNOWN
                << " nMIXED:" << nMIXED
                << " nINSIDE:" << nINSIDE
                << " nOUTSIDE:" << nOUTSIDE
                << endl;
        }
    }

    return getVolumeType(0, sample);
}


template<class Type>
template<class CompareOp>
void Foam::dynamicIndexedOctree<Type>::findNear
(
    const scalar nearDist,
    const dynamicIndexedOctree<Type>& tree2,
    CompareOp& cop
) const
{
    findNear
    (
        nearDist,
        true,
        *this,
        nodePlusOctant(0, 0),
        bb(),
        tree2,
        nodePlusOctant(0, 0),
        tree2.bb(),
        cop
    );
}


template<class Type>
bool Foam::dynamicIndexedOctree<Type>::insert(label startIndex, label endIndex)
{
    if (startIndex == endIndex)
    {
        return false;
    }

    if (nodes_.empty())
    {
        contents_.append
        (
            autoPtr<DynamicList<label> >
            (
                new DynamicList<label>(1)
            )
        );

        contents_[0]().append(0);

        // Create topnode.
        node topNode = divide(bb_, 0, -1, 0);

        nodes_.append(topNode);

        startIndex++;
    }

    bool success = true;

    for (label pI = startIndex; pI < endIndex; ++pI)
    {
        label nLevels = 1;

        if (!insertIndex(0, pI, nLevels))
        {
            success = false;
        }

        nLevelsMax_ = max(nLevels, nLevelsMax_);
    }

    return success;
}


template<class Type>
bool Foam::dynamicIndexedOctree<Type>::insertIndex
(
    const label nodIndex,
    const label index,
    label& nLevels
)
{
    bool shapeInserted = false;

    for (direction octant = 0; octant < 8; octant++)
    {
        const labelBits& subNodeLabel = nodes_[nodIndex].subNodes_[octant];

        if (isNode(subNodeLabel))
        {
            const treeBoundBox& subBb = nodes_[getNode(subNodeLabel)].bb_;

            if (shapes().overlaps(index, subBb))
            {
                nLevels++;

                if (insertIndex(getNode(subNodeLabel), index, nLevels))
                {
                    shapeInserted = true;
                }
            }
        }
        else if (isContent(subNodeLabel))
        {
            const treeBoundBox subBb = nodes_[nodIndex].bb_.subBbox(octant);

            if (shapes().overlaps(index, subBb))
            {
                const label contentI = getContent(subNodeLabel);

                contents_[contentI]().append(index);

                recursiveSubDivision
                (
                    subBb,
                    contentI,
                    nodIndex,
                    octant,
                    nLevels
                );

                shapeInserted = true;
            }
        }
        else
        {
            const treeBoundBox subBb = nodes_[nodIndex].bb_.subBbox(octant);

            if (shapes().overlaps(index, subBb))
            {
                label sz = contents_.size();

                contents_.append
                (
                    autoPtr<DynamicList<label> >(new DynamicList<label>(1))
                );

                contents_[sz]().append(index);

                nodes_[nodIndex].subNodes_[octant]
                    = contentPlusOctant(sz, octant);
            }

            shapeInserted = true;
        }
    }

    return shapeInserted;
}


template<class Type>
bool Foam::dynamicIndexedOctree<Type>::remove(const label index)
{
    if (nodes_.empty())
    {
        return false;
    }

    removeIndex(0, index);

    return true;
}


template<class Type>
Foam::label Foam::dynamicIndexedOctree<Type>::removeIndex
(
    const label nodIndex,
    const label index
)
{
    label totalContents = 0;

    for (direction octant = 0; octant < 8; octant++)
    {
        const labelBits& subNodeLabel = nodes_[nodIndex].subNodes_[octant];

        if (isNode(subNodeLabel))
        {
            const treeBoundBox& subBb = nodes_[getNode(subNodeLabel)].bb_;

            if (shapes().overlaps(index, subBb))
            {
                // Recursive function.
                label childContentsSize
                    = removeIndex(getNode(subNodeLabel), index);

                totalContents += childContentsSize;

                if (childContentsSize == 0)
                {
                    nodes_[nodIndex].subNodes_[octant]
                        = emptyPlusOctant(octant);
                }
            }
            else
            {
                // Increment this so that the node will not be marked as empty
                totalContents++;
            }
        }
        else if (isContent(subNodeLabel))
        {
            const treeBoundBox subBb = nodes_[nodIndex].bb_.subBbox(octant);

            const label contentI = getContent(subNodeLabel);

            if (shapes().overlaps(index, subBb))
            {
                DynamicList<label>& contentList = contents_[contentI]();

                DynamicList<label> newContent(contentList.size());

                forAll(contentList, pI)
                {
                    const label oldIndex = contentList[pI];

                    if (oldIndex != index)
                    {
                        newContent.append(oldIndex);
                    }
                }

                newContent.shrink();

                if (newContent.size() == 0)
                {
                    // Set to empty.
                    nodes_[nodIndex].subNodes_[octant]
                        = emptyPlusOctant(octant);
                }

                contentList.transfer(newContent);
            }

            totalContents += contents_[contentI]().size();
        }
        else
        {
            // Empty, do nothing.
        }
    }

    return totalContents;
}


// Print contents of nodeI
template<class Type>
void Foam::dynamicIndexedOctree<Type>::print
(
    prefixOSstream& os,
    const bool printContents,
    const label nodeI
) const
{
    const node& nod = nodes_[nodeI];
    const treeBoundBox& bb = nod.bb_;

    os  << "nodeI:" << nodeI << " bb:" << bb << nl
        << "parent:" << nod.parent_ << nl
        << "n:" << countElements(nodePlusOctant(nodeI, 0)) << nl;

    for (direction octant = 0; octant < nod.subNodes_.size(); octant++)
    {
        const treeBoundBox subBb(bb.subBbox(octant));

        labelBits index = nod.subNodes_[octant];

        if (isNode(index))
        {
            // tree node.
            label subNodeI = getNode(index);

            os  << "octant:" << octant
                << " node: n:" << countElements(index)
                << " bb:" << subBb << endl;

            string oldPrefix = os.prefix();
            os.prefix() = "  " + oldPrefix;

            print(os, printContents, subNodeI);

            os.prefix() = oldPrefix;
        }
        else if (isContent(index))
        {
            const labelList& indices = contents_[getContent(index)];

            if (false) //debug)
            {
                writeOBJ(nodeI, octant);
            }

            os  << "octant:" << octant
                << " content: n:" << indices.size()
                << " bb:" << subBb;

            if (printContents)
            {
                os << " contents:";
                forAll(indices, j)
                {
                    os  << ' ' << indices[j];
                }
            }
            os  << endl;
        }
        else
        {
            if (false) //debug)
            {
                writeOBJ(nodeI, octant);
            }

            os  << "octant:" << octant << " empty:" << subBb << endl;
        }
    }
}


template<class Type>
void Foam::dynamicIndexedOctree<Type>::writeTreeInfo() const
{
    label nEntries = 0;
    forAll(contents_, i)
    {
        nEntries += contents_[i]().size();
    }

    Pout<< "indexedOctree<Type>::indexedOctree"
        << " : finished construction of tree of:" << shapes().typeName
        << nl
        << "    bounding box:     " << this->bb() << nl
        << "    shapes:           " << shapes().size() << nl
        << "    treeNodes:        " << nodes_.size() << nl
        << "    nEntries:         " << nEntries << nl
        << "    levels/maxLevels: " << nLevelsMax_ << "/" << maxLevels_ << nl
        << "    minSize:          " << minSize_ << nl
        << "        per treeLeaf:         "
        << scalar(nEntries)/contents_.size() << nl
        << "        per shape (duplicity):"
        << scalar(nEntries)/shapes().size() << nl
        << endl;
}


// Print contents of nodeI
template<class Type>
bool Foam::dynamicIndexedOctree<Type>::write(Ostream& os) const
{
    os << *this;

    return os.good();
}


template<class Type>
Foam::Ostream&
Foam::operator<<(Ostream& os, const dynamicIndexedOctree<Type>& t)
{
    os  << t.bb() << token::SPACE << t.nodes() << endl;

    forAll(t.contents(), cI)
    {
        os << t.contents()[cI]() << endl;
    }

    return os;
}


// ************************************************************************* //
