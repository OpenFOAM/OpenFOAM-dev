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

Class
    Foam::surfaceSets

Description
    Various utilities to handle sets relating mesh to surface.
    Note: work in progress. Used in meshing tools.

SourceFiles
    surfaceSets.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSets_H
#define surfaceSets_H

#include "label.H"
#include "scalar.H"
#include "fileName.H"
#include "pointField.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class triSurface;
class triSurfaceSearch;
class pointSet;

/*---------------------------------------------------------------------------*\
                           Class surfaceSets Declaration
\*---------------------------------------------------------------------------*/

class surfaceSets
{
    //- Visibility of points: from inside cells, from outside cells or
    //  from both
    enum pointStatus
    {
        NOTSET,
        INSIDE,
        MIXED,
        OUTSIDE
    };


    // Static

        ////- Min length of connected mesh edges
        // static scalar minEdgeLen
        //(
        //    const primitiveMesh& mesh,
        //    const label pointi
        //);
        //
        ////- Returns true if cell uses at least one selected point
        // static bool usesPoint
        //(
        //    const primitiveMesh& mesh,
        //    const boolList& selectedPoint,
        //    const label celli
        //);

        ////- Find points used only by internalCells and remove from
        ////  internalCells all
        ////  cells which don't use any internal point. Return number of
        ////  removed cells.
        // static label removeHangingCells
        //(
        //    const primitiveMesh&,
        //    const triSurfaceSearch& querySurf,
        //    labelHashSet& internalCells
        //);

        ////- Select all points out of pointSet where the distance to
        ////  the surface is less than a factor times a local length
        ////  scale (minimum length of connected edges)
        // static void getNearPoints
        //(
        //    const primitiveMesh& mesh,
        //    const triSurface& surf,
        //    const triSurfaceSearch& querySurf,
        //    const scalar edgeFactor,
        //    const pointSet& candidateSet,
        //    pointSet& nearPointSet
        //);

public:

    // Static Functions

        //- Divide cells into cut,inside and outside
        //  nCutLayers>0 : remove cutCells (set to type inside) if further
        //  than nCutLayers away from outside type cell.
        static void getSurfaceSets
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const triSurface& surf,
            const triSurfaceSearch& querySurf,
            const pointField& outsidePts,

            const label nCutLayers,

            labelHashSet& inside,
            labelHashSet& outside,
            labelHashSet& cut
        );

        //- Get cells using points on 'outside' only
        static labelHashSet getHangingCells
        (
            const primitiveMesh& mesh,
            const labelHashSet& internalCells
        );


//        //- Write cell sets with cells 'inside' and 'outside' surface.
//        static void writeSurfaceSets
//        (
//            const polyMesh& mesh,
//            const fileName& surfName,
//            const triSurface& surf,
//            const triSurfaceSearch& querySurf,
//            const pointField& outsidePts,
//            const scalar edgeFactor
//        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
