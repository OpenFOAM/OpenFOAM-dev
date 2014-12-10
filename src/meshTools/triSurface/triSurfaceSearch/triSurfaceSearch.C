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

#include "triSurfaceSearch.H"
#include "triSurface.H"
#include "PatchTools.H"
#include "volumeType.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triSurfaceSearch::checkUniqueHit
(
    const pointIndexHit& currHit,
    const DynamicList<pointIndexHit, 1, 1>& hits,
    const vector& lineVec
) const
{
    // Classify the hit
    label nearType = -1;
    label nearLabel = -1;

    const labelledTri& f = surface()[currHit.index()];

    f.nearestPointClassify
    (
        currHit.hitPoint(),
        surface().points(),
        nearType,
        nearLabel
    );

    if (nearType == triPointRef::POINT)
    {
        // near point

        const label nearPointI = f[nearLabel];

        const labelList& pointFaces =
            surface().pointFaces()[surface().meshPointMap()[nearPointI]];

        forAll(pointFaces, pI)
        {
            const label pointFaceI = pointFaces[pI];

            if (pointFaceI != currHit.index())
            {
                forAll(hits, hI)
                {
                    const pointIndexHit& hit = hits[hI];

                    if (hit.index() == pointFaceI)
                    {
                        return false;
                    }
                }
            }
        }
    }
    else if (nearType == triPointRef::EDGE)
    {
        // near edge
        // check if the other face of the edge is already hit

        const labelList& fEdges = surface().faceEdges()[currHit.index()];

        const label edgeI = fEdges[nearLabel];

        const labelList& edgeFaces = surface().edgeFaces()[edgeI];

        forAll(edgeFaces, fI)
        {
            const label edgeFaceI = edgeFaces[fI];

            if (edgeFaceI != currHit.index())
            {
                forAll(hits, hI)
                {
                    const pointIndexHit& hit = hits[hI];

                    if (hit.index() == edgeFaceI)
                    {
                        // Check normals
                        const vector currHitNormal =
                            surface().faceNormals()[currHit.index()];

                        const vector existingHitNormal =
                            surface().faceNormals()[edgeFaceI];

                        const label signCurrHit =
                            pos(currHitNormal & lineVec);

                        const label signExistingHit =
                            pos(existingHitNormal & lineVec);

                        if (signCurrHit == signExistingHit)
                        {
                            return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceSearch::triSurfaceSearch(const triSurface& surface)
:
    surface_(surface),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    maxTreeDepth_(10),
    treePtr_(NULL)
{}


Foam::triSurfaceSearch::triSurfaceSearch
(
    const triSurface& surface,
    const dictionary& dict
)
:
    surface_(surface),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    maxTreeDepth_(10),
    treePtr_(NULL)
{
    // Have optional non-standard search tolerance for gappy surfaces.
    if (dict.readIfPresent("tolerance", tolerance_) && tolerance_ > 0)
    {
        Info<< "    using intersection tolerance " << tolerance_ << endl;
    }

    // Have optional non-standard tree-depth to limit storage.
    if (dict.readIfPresent("maxTreeDepth", maxTreeDepth_) && maxTreeDepth_ > 0)
    {
        Info<< "    using maximum tree depth " << maxTreeDepth_ << endl;
    }
}


Foam::triSurfaceSearch::triSurfaceSearch
(
    const triSurface& surface,
    const scalar tolerance,
    const label maxTreeDepth
)
:
    surface_(surface),
    tolerance_(tolerance),
    maxTreeDepth_(maxTreeDepth),
    treePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceSearch::~triSurfaceSearch()
{
    clearOut();
}


void Foam::triSurfaceSearch::clearOut()
{
    treePtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataTriSurface>&
Foam::triSurfaceSearch::tree() const
{
    if (treePtr_.empty())
    {
        // Calculate bb without constructing local point numbering.
        treeBoundBox bb(vector::zero, vector::zero);

        if (surface().size())
        {
            label nPoints;
            PatchTools::calcBounds(surface(), bb, nPoints);

            if (nPoints != surface().points().size())
            {
                WarningIn("triSurfaceSearch::tree() const")
                    << "Surface does not have compact point numbering."
                    << " Of " << surface().points().size()
                    << " only " << nPoints
                    << " are used. This might give problems in some routines."
                    << endl;
            }

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
        indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

        treePtr_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(true, surface_, tolerance_),
                bb,
                maxTreeDepth_,  // maxLevel
                10,             // leafsize
                3.0             // duplicity
            )
        );

        indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
    }

    return treePtr_();
}


// Determine inside/outside for samples
Foam::boolList Foam::triSurfaceSearch::calcInside
(
    const pointField& samples
) const
{
    boolList inside(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        if (!tree().bb().contains(sample))
        {
            inside[sampleI] = false;
        }
        else if (tree().getVolumeType(sample) == volumeType::INSIDE)
        {
            inside[sampleI] = true;
        }
        else
        {
            inside[sampleI] = false;
        }
    }
    return inside;
}


void Foam::triSurfaceSearch::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(samples.size());

    forAll(samples, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findNearest
        (
            samples[i],
            nearestDistSqr[i],
            treeDataTriSurface::findNearestOp(octree)
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


Foam::pointIndexHit Foam::triSurfaceSearch::nearest
(
    const point& pt,
    const vector& span
)
const
{
    const scalar nearestDistSqr = 0.25*magSqr(span);

    return tree().findNearest(pt, nearestDistSqr);
}


void Foam::triSurfaceSearch::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    forAll(start, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findLine
        (
            start[i],
            end[i]
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceSearch::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    forAll(start, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findLineAny
        (
            start[i],
            end[i]
        );
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceSearch::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    // Work array
    DynamicList<pointIndexHit, 1, 1> hits;

    DynamicList<label> shapeMask;

    treeDataTriSurface::findAllIntersectOp allIntersectOp(octree, shapeMask);

    forAll(start, pointI)
    {
        hits.clear();
        shapeMask.clear();

        while (true)
        {
            // See if any intersection between pt and end
            pointIndexHit inter = octree.findLine
            (
                start[pointI],
                end[pointI],
                allIntersectOp
            );

            if (inter.hit())
            {
                vector lineVec = end[pointI] - start[pointI];
                lineVec /= mag(lineVec) + VSMALL;

                if
                (
                    checkUniqueHit
                    (
                        inter,
                        hits,
                        lineVec
                    )
                )
                {
                    hits.append(inter);
                }

                shapeMask.append(inter.index());
            }
            else
            {
                break;
            }
        }

        info[pointI].transfer(hits);
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


// ************************************************************************* //
