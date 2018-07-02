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

#include "triSurfaceRegionSearch.H"
#include "indexedOctree.H"
#include "triSurface.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceRegionSearch::triSurfaceRegionSearch(const triSurface& surface)
:
    triSurfaceSearch(surface),
    indirectRegionPatches_(),
    treeByRegion_()
{}


Foam::triSurfaceRegionSearch::triSurfaceRegionSearch
(
    const triSurface& surface,
    const dictionary& dict
)
:
    triSurfaceSearch(surface, dict),
    indirectRegionPatches_(),
    treeByRegion_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceRegionSearch::~triSurfaceRegionSearch()
{
    clearOut();
}


void Foam::triSurfaceRegionSearch::clearOut()
{
    triSurfaceSearch::clearOut();
    treeByRegion_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::PtrList<Foam::triSurfaceRegionSearch::treeType>&
Foam::triSurfaceRegionSearch::treeByRegion() const
{
    if (treeByRegion_.empty())
    {
        Map<label> regionSizes;
        forAll(surface(), fI)
        {
            const label regionI = surface()[fI].region();

            regionSizes(regionI)++;
        }

        label nRegions = regionSizes.size();

        indirectRegionPatches_.setSize(nRegions);
        treeByRegion_.setSize(nRegions);

        labelListList regionsAddressing(nRegions);

        forAll(regionsAddressing, regionI)
        {
            regionsAddressing[regionI] = labelList(regionSizes[regionI], -1);
        }

        labelList nFacesInRegions(nRegions, 0);

        forAll(surface(), fI)
        {
            const label regionI = surface()[fI].region();

            regionsAddressing[regionI][nFacesInRegions[regionI]++] = fI;
        }

        forAll(regionsAddressing, regionI)
        {
            scalar oldTol = treeType::perturbTol();
            treeType::perturbTol() = tolerance();

            indirectRegionPatches_.set
            (
                regionI,
                new indirectTriSurface
                (
                    IndirectList<labelledTri>
                    (
                        surface(),
                        regionsAddressing[regionI]
                    ),
                    surface().points()
                )
            );

            // Calculate bb without constructing local point numbering.
            treeBoundBox bb(Zero, Zero);

            if (indirectRegionPatches_[regionI].size())
            {
                label nPoints;
                PatchTools::calcBounds
                (
                    indirectRegionPatches_[regionI],
                    bb,
                    nPoints
                );

    //            if (nPoints != surface().points().size())
    //            {
    //                WarningInFunction
    //                    << "Surface does not have compact point numbering. "
    //                    << "Of " << surface().points().size()
    //                    << " only " << nPoints
    //                    << " are used."
    //                    << " This might give problems in some routines."
    //                    << endl;
    //            }

                // Slightly extended bb. Slightly off-centred just so
                // on symmetric geometry there are fewer face/edge
                // aligned items.
                bb = bb.extend(1e-4);
            }

            treeByRegion_.set
            (
                regionI,
                new treeType
                (
                    treeDataIndirectTriSurface
                    (
                        true,
                        indirectRegionPatches_[regionI],
                        tolerance()
                    ),
                    bb,
                    maxTreeDepth(),  // maxLevel
                    10,              // leafsize
                    3.0              // duplicity
                )
            );

            treeType::perturbTol() = oldTol;
        }
    }

    return treeByRegion_;
}


void Foam::triSurfaceRegionSearch::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    if (regionIndices.empty())
    {
        triSurfaceSearch::findNearest(samples, nearestDistSqr, info);
    }
    else
    {
        scalar oldTol = treeType::perturbTol();
        treeType::perturbTol() = tolerance();

        const PtrList<treeType>& octrees = treeByRegion();

        info.setSize(samples.size());

        forAll(octrees, treeI)
        {
            if (findIndex(regionIndices, treeI) == -1)
            {
                continue;
            }

            const treeType& octree = octrees[treeI];

            forAll(samples, i)
            {
//                if (!octree.bb().contains(samples[i]))
//                {
//                    continue;
//                }

                pointIndexHit currentRegionHit = octree.findNearest
                (
                    samples[i],
                    nearestDistSqr[i],
                    treeDataIndirectTriSurface::findNearestOp(octree)
                );

                if
                (
                    currentRegionHit.hit()
                 &&
                    (
                        !info[i].hit()
                     ||
                        (
                            magSqr(currentRegionHit.hitPoint() - samples[i])
                          < magSqr(info[i].hitPoint() - samples[i])
                        )
                    )
                )
                {
                    info[i] = currentRegionHit;
                }
            }
        }

        treeType::perturbTol() = oldTol;
    }
}


// ************************************************************************* //
