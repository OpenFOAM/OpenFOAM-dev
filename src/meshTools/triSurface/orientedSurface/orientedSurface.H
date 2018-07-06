/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::orientedSurface

Description
    Given point flip all faces such that normals point in same direction.

SourceFiles
    orientedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef orientedSurface_H
#define orientedSurface_H

#include "triSurface.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurfaceSearch;

/*---------------------------------------------------------------------------*\
                       Class orientedSurface Declaration
\*---------------------------------------------------------------------------*/

class orientedSurface
:
    public triSurface
{
public:

    // Data types

        //- Enumeration listing whether face needs to be flipped.
        enum sideStat
        {
            UNVISITED,
            FLIP,
            NOFLIP
        };

private:

    // Private Member Functions

        //- Return true if edge is used in opposite order in faces
        static bool consistentEdge
        (
            const edge& e,
            const triSurface::FaceType& f0,
            const triSurface::FaceType& f1
        );

        //- From changed faces get the changed edges
        static labelList faceToEdge
        (
            const triSurface&,
            const labelList& changedFaces
        );

        //- From changed edges check the orientation of the connected faces
        //  and flip them. Return changed faces.
        static labelList edgeToFace
        (
            const triSurface&,
            const labelList& changedEdges,
            labelList& flip
        );

        //- Walk from face across connected faces. Change orientation to be
        //  consistent with startFacei.
        static void walkSurface
        (
            const triSurface& s,
            const label startFacei,
            labelList& flipState
        );

        //- Given nearest point and face check orientation to nearest face
        //  and flip if necessary (only marked in flipState) and propagate.
        static void propagateOrientation
        (
            const triSurface&,
            const point& outsidePoint,
            const bool orientOutside,
            const label nearestFacei,
            const point& nearestPt,
            labelList& flipState
        );

        //- Find a face on zoneI and count number of intersections to determine
        //  orientation
        static void findZoneSide
        (
            const triSurfaceSearch& surfSearches,
            const labelList& faceZone,
            const label zoneI,
            const point& visiblePoint,
            label& zoneFacei,
            bool& isOutside
        );

        //- Given flipState reverse triangles of *this. Return true if
        //  anything flipped.
        static bool flipSurface(triSurface& s, const labelList& flipState);

        //- Make surface surface has consistent orientation across connected
        //  triangles.
        static bool orientConsistent(triSurface& s);


public:

    ClassName("orientedSurface");


    // Constructors

        //- Construct null
        orientedSurface();

        //- Construct from triSurface and sample point which is either
        //  outside (orientOutside = true) or inside (orientOutside = false).
        //  Uses orient.
        orientedSurface
        (
            const triSurface&,
            const point& samplePoint,
            const bool orientOutside = true
        );

        //- Construct from triSurface. Calculates outside point as being
        //  outside the bounding box of the surface.
        orientedSurface(const triSurface&, const bool orientOutside = true);


    // Member Functions

        //- Flip faces such that normals are consistent with point:
        //  orientOutside=true : point outside surface
        //  orientOutside=false : point inside surface
        //  Bases orientation on normal on nearest point (linear search) and
        //  walks to rest. Surface needs to be manifold.
        static bool orient(triSurface&, const point&, const bool orientOutside);

        //- Flip faces such that normals are consistent with point:
        //  orientOutside=true : point outside surface
        //  orientOutside=false : point inside surface
        //  Uses intersection count to orient. Handles open surfaces.
        static bool orient
        (
            triSurface& s,
            const triSurfaceSearch& querySurf,
            const point& samplePoint,
            const bool orientOutside
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
