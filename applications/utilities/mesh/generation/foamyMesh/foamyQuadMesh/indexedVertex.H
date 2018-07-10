/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    CGAL::indexedVertex

Description
    An indexed form of CGAL::Triangulation_vertex_base_2<K> used to keep
    track of the vertices in the triangulation.

\*---------------------------------------------------------------------------*/

#ifndef indexedVertex_H
#define indexedVertex_H

#include "CGAL/Triangulation_2.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{

// Forward declaration of friend functions and operators

template<class Gt, class Vb>
class indexedVertex;

template<class Gt, class Vb>
bool pointPair
(
    const indexedVertex<Gt, Vb>& v0,
    const indexedVertex<Gt, Vb>& v1
);

template<class Gt, class Vb>
bool boundaryTriangle
(
    const indexedVertex<Gt, Vb>& v0,
    const indexedVertex<Gt, Vb>& v1,
    const indexedVertex<Gt, Vb>& v2
);

template<class Gt, class Vb>
bool outsideTriangle
(
    const indexedVertex<Gt, Vb>& v0,
    const indexedVertex<Gt, Vb>& v1,
    const indexedVertex<Gt, Vb>& v2
);

/*---------------------------------------------------------------------------*\
                       Class indexedVertex Declaration
\*---------------------------------------------------------------------------*/

template<class Gt, class Vb=CGAL::Triangulation_vertex_base_2<Gt>>
class indexedVertex
:
    public Vb
{
    // Private data

        //- The index for this triangle vertex
        int index_;

        //- Index of pair-point :
        //  NEAR_BOUNDARY_POINT : internal near boundary point.
        //  INTERNAL_POINT      : internal point.
        //  FAR_POINT           : far-point.
        //  >= 0                : part of point-pair. Index of other point.
        //                        Lowest numbered is inside one (master).
        int type_;


public:

    enum pointTypes
    {
        NEAR_BOUNDARY_POINT = -4,
        INTERNAL_POINT      = -3,
        MIRROR_POINT        = -2,
        FAR_POINT           = -1
    };

    typedef typename Vb::Vertex_handle      Vertex_handle;
    typedef typename Vb::Face_handle        Face_handle;
    typedef typename Vb::Point              Point;

    template<class TDS2>
    struct Rebind_TDS
    {
        typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
        typedef indexedVertex<Gt,Vb2>                           Other;
    };


    // Constructors

        inline indexedVertex();

        inline indexedVertex(const Point& p);

        inline indexedVertex(const Point& p, const int index, const int& type);

        inline indexedVertex(const Point& p, Face_handle f);

        inline indexedVertex(Face_handle f);


    // Member Functions

        inline int& index();

        inline int index() const;

        inline int& type();

        inline int type() const;

        //- Is point a far-point
        inline bool farPoint() const;

        //- Is point internal, i.e. not on boundary
        inline bool internalPoint() const;

        //- Is point internal and near the boundary
        inline bool nearBoundary() const;

        //- Set the point to be near the boundary
        inline void setNearBoundary();

        //- Is point a mirror point
        inline bool mirrorPoint() const;

        //- Either master or slave of pointPair.
        inline bool pairPoint() const;

        //- Master of a pointPair is the lowest numbered one.
        inline bool ppMaster() const;

        //- Slave of a pointPair is the highest numbered one.
        inline bool ppSlave() const;

        //- Either original internal point or master of pointPair.
        inline bool internalOrBoundaryPoint() const;

        //- Is point near the boundary or part of the boundary definition
        inline bool nearOrOnBoundary() const;


    // Friend Functions

        //- Do the two given vertices consitute a boundary point-pair
        friend bool pointPair <Gt, Vb>
        (
            const indexedVertex<Gt, Vb>& v0,
            const indexedVertex<Gt, Vb>& v1
        );

        //- Do the three given vertices consitute a boundary triangle
        friend bool boundaryTriangle <Gt, Vb>
        (
            const indexedVertex<Gt, Vb>& v0,
            const indexedVertex<Gt, Vb>& v1,
            const indexedVertex<Gt, Vb>& v2
        );

        //- Do the three given vertices consitute an outside triangle
        friend bool outsideTriangle <Gt, Vb>
        (
            const indexedVertex<Gt, Vb>& v0,
            const indexedVertex<Gt, Vb>& v1,
            const indexedVertex<Gt, Vb>& v2
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "indexedVertexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
