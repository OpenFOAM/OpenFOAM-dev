/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Gt, class Fb>
inline CGAL::indexedFace<Gt, Fb>::indexedFace()
:
    Fb(),
    index_(CHANGED)
{}


template<class Gt, class Fb>
inline CGAL::indexedFace<Gt, Fb>::indexedFace
(
    Vertex_handle v0,
    Vertex_handle v1,
    Vertex_handle v2
)
:
    Fb(v0, v1, v2),
    index_(CHANGED)
{}


template<class Gt, class Fb>
inline CGAL::indexedFace<Gt, Fb>::indexedFace
(
    Vertex_handle v0,
    Vertex_handle v1,
    Vertex_handle v2,
    Face_handle n0,
    Face_handle n1,
    Face_handle n2
)
:
    Fb(v0, v1, v2, n0, n1, n2),
    index_(CHANGED)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Gt, class Fb>
inline void CGAL::indexedFace<Gt, Fb>::set_vertex(int i, Vertex_handle v)
{
    index_ = CHANGED;
    Fb::set_vertex(i, v);
}


template<class Gt, class Fb>
inline void CGAL::indexedFace<Gt, Fb>::set_vertices()
{
    index_ = CHANGED;
    Fb::set_vertices();
}


template<class Gt, class Fb>
inline void CGAL::indexedFace<Gt, Fb>::set_vertices
(
    Vertex_handle v0,
    Vertex_handle v1,
    Vertex_handle v2
)
{
    index_ = CHANGED;
    Fb::set_vertices(v0, v1, v2);
}


template<class Gt, class Fb>
inline int& CGAL::indexedFace<Gt, Fb>::faceIndex()
{
    return index_;
}


template<class Gt, class Fb>
inline int CGAL::indexedFace<Gt, Fb>::faceIndex() const
{
    return index_;
}


// ************************************************************************* //
