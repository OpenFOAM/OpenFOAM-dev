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

Class
    Foam::triSurfaceRegionSearch

Description
    Helper class to search on triSurface. Creates an octree for each region of
    the surface and only searches on the specified regions.

SourceFiles
    triSurfaceRegionSearch.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceRegionSearch_H
#define triSurfaceRegionSearch_H

#include "pointField.H"
#include "pointIndexHit.H"
#include "triSurfaceSearch.H"
#include "labelledTri.H"
#include "IndirectList.H"
#include "PtrList.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class triSurfaceRegionSearch Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceRegionSearch
:
    public triSurfaceSearch
{
    // Private typedefs

        typedef PrimitivePatch
        <
            labelledTri,
            IndirectList,
            const pointField&
        > indirectTriSurface;

        typedef treeDataPrimitivePatch<indirectTriSurface>
            treeDataIndirectTriSurface;

        typedef indexedOctree<treeDataIndirectTriSurface> treeType;


    // Private data

        //- Surface is split into patches by region
        mutable PtrList<indirectTriSurface> indirectRegionPatches_;

        //- Search tree for each region
        mutable PtrList<treeType> treeByRegion_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        triSurfaceRegionSearch(const triSurfaceRegionSearch&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceRegionSearch&);


public:

    // Constructors

        //- Construct from surface. Holds reference to surface!
        explicit triSurfaceRegionSearch(const triSurface&);

        //- Construct from surface and dictionary. Holds reference to surface!
        triSurfaceRegionSearch(const triSurface&, const dictionary& dict);


    //- Destructor
    ~triSurfaceRegionSearch();

        //- Clear storage
        void clearOut();


    // Member Functions

        // Access

            //- Demand driven construction of octree for each region.
            //  @todo Currently creates a tree for each region; could optimise
            //        by only constructing trees when they are in regionIndices
            const PtrList<treeType>& treeByRegion() const;


        // Query

            //- Find the nearest point on the surface out of the regions
            //  supplied in the list regionIndices. Ignores regions that are
            //  not specified
            void findNearest
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                List<pointIndexHit>& info
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
