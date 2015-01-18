/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    indexedFace

Description
    An indexed form of CGAL::Triangulation_face_base_2<K> used to keep
    track of the vertices in the triangulation.

\*---------------------------------------------------------------------------*/

#ifndef indexedFace_H
#define indexedFace_H

#include "CGAL/Triangulation_2.h"
#include "indexedVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{

/*---------------------------------------------------------------------------*\
                       Class indexedFace Declaration
\*---------------------------------------------------------------------------*/

template<class Gt, class Fb=CGAL::Triangulation_face_base_2<Gt> >
class indexedFace
:
    public Fb
{
    // Private data

        //- The index for this triangle face
        //  -1: triangle and changed and associated data needs updating
        // >=0: index of triangles, set by external update algorithm
        int index_;


public:

    enum faceTypes
    {
        UNCHANGED           = 0,
        CHANGED             = -1,
        SAVE_CHANGED        = -2
    };

    typedef typename Fb::Vertex_handle                   Vertex_handle;
    typedef typename Fb::Face_handle                     Face_handle;

    template<typename TDS2>
    struct Rebind_TDS
    {
        typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef indexedFace<Gt, Fb2>  Other;
    };


    // Constructors

        inline indexedFace();

        inline indexedFace
        (
            Vertex_handle v0,
            Vertex_handle v1,
            Vertex_handle v2
        );

        inline indexedFace
        (
            Vertex_handle v0,
            Vertex_handle v1,
            Vertex_handle v2,
            Face_handle n0,
            Face_handle n1,
            Face_handle n2
        );


    // Member Functions

        inline void set_vertex(int i, Vertex_handle v);

        inline void set_vertices();

        inline void set_vertices
        (
            Vertex_handle v0,
            Vertex_handle v1,
            Vertex_handle v2
        );

        inline int& faceIndex();

        inline int faceIndex() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "indexedFaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
