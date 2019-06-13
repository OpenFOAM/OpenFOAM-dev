/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::edgeSurface

Description
    Description of surface in form of 'cloud of edges'.

    The 'cloud of edges':
    - points
    - edges
    - faceEdges
    - parentEdge (edge on surface this edge originates from)
    and nothing more.

    (pointEdges constructed from above data)

    Constructed from triSurface and surfaceIntersection. (uses localPoints
    of surface of course)

    Used to easily insert cuts and split faces.

Note
    - points with surface (local)points first, intersection points last
    - edges with (split) surface edges first, intersection edges last.

SourceFiles
    edgeSurface.C

\*---------------------------------------------------------------------------*/

#ifndef edgeSurface_H
#define edgeSurface_H

#include "edgeList.H"
#include "labelList.H"
#include "pointField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;
class surfaceIntersection;

/*---------------------------------------------------------------------------*\
                           Class edgeSurface Declaration
\*---------------------------------------------------------------------------*/

class edgeSurface
{
    // Private Data

        //- All points (0 .. nSurfacePoints_-1 are points from surface)
        pointField points_;

        label nSurfacePoints_;

        //- All edges (0 .. nSurfaceEdges_-1 are (possibly split) surface edges)
        edgeList edges_;

        label nSurfaceEdges_;

        //- Original surface edge. Valid only surfaceEdges.
        labelList parentEdges_;

        //- From face to our edges_
        labelListList faceEdges_;


        //- Constructed from above: pointEdges
        labelListList pointEdges_;


    // Private Member Functions

        //- Dump edges in obj format
        static void writeOBJ(const pointField&, const edgeList&, Ostream&);

        //- Dump selected edges in obj format
        static void writeOBJ
        (
            const pointField&,
            const edgeList&,
            const labelList&,
            Ostream&
        );

        //- Calculate pointEdges
        void calcPointEdges();


public:

    ClassName("edgeSurface");

    // Constructors

        //- Construct from surface and intersection description
        edgeSurface
        (
            const triSurface& surf,
            const bool isFirstSurface,
            const surfaceIntersection& inter
        );


    // Member Functions

        // Access

            const pointField& points() const
            {
                return points_;
            }

            label nSurfacePoints() const
            {
                return nSurfacePoints_;
            }

            const edgeList& edges() const
            {
                return edges_;
            }

            label nSurfaceEdges() const
            {
                return nSurfaceEdges_;
            }

            bool isSurfaceEdge(const label edgeI) const
            {
                return edgeI < nSurfaceEdges_;
            }

            //- Parent edge (original surface edge this edge came from).
            //  Valid only for edgeI < nSurfaceEdges_.
            label parentEdge(const label edgeI) const
            {
                if (edgeI < nSurfaceEdges_)
                {
                    return parentEdges_[edgeI];
                }
                else
                {
                    FatalErrorInFunction
                        << "Trying to get parent (i.e. surface) edge for"
                        << " intersection edge " << edgeI
                        << abort(FatalError);
                    return -1;
                }
            }

            //- From face to our edges_
            const labelListList& faceEdges() const
            {
                return faceEdges_;
            }

            //- Point to edge addressing
            const labelListList& pointEdges() const
            {
                return pointEdges_;
            }


        // Edit

            //- Add intersection edges to a face. Used for connecting
            //  floating intersection on face to rest of face.
            void addIntersectionEdges(const label facei, const edgeList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
