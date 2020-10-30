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

template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex()
:
    Vb(),
    index_(INTERNAL_POINT),
    type_(INTERNAL_POINT)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(const Point& p)
:
    Vb(p),
    index_(INTERNAL_POINT),
    type_(INTERNAL_POINT)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex
(
    const Point& p,
    const int index,
    const int& type
)
:
    Vb(p),
    index_(index),
    type_(type)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(const Point& p, Face_handle f)
:
    Vb(f, p),
    index_(INTERNAL_POINT),
    type_(INTERNAL_POINT)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(Face_handle f)
:
    Vb(f),
    index_(INTERNAL_POINT),
    type_(INTERNAL_POINT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Gt, class Vb>
inline int& CGAL::indexedVertex<Gt, Vb>::index()
{
    return index_;
}


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::index() const
{
    return index_;
}


template<class Gt, class Vb>
inline int& CGAL::indexedVertex<Gt, Vb>::type()
{
    return type_;
}


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::type() const
{
    return type_;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::farPoint() const
{
    return type_ == FAR_POINT;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalPoint() const
{
    return type_ <= INTERNAL_POINT;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::nearBoundary() const
{
    return type_ == NEAR_BOUNDARY_POINT;
}


template<class Gt, class Vb>
inline void CGAL::indexedVertex<Gt, Vb>::setNearBoundary()
{
    type_ = NEAR_BOUNDARY_POINT;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::mirrorPoint() const
{
    return type_ == MIRROR_POINT;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::pairPoint() const
{
    return type_ >= 0;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::ppMaster() const
{
    if (type_ > index_)
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::ppSlave() const
{
    if (type_ >= 0 && type_ < index_)
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalOrBoundaryPoint() const
{
    return internalPoint() || ppMaster();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::nearOrOnBoundary() const
{
    return pairPoint() || mirrorPoint() || nearBoundary();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Gt, class Vb>
bool CGAL::pointPair
(
    const indexedVertex<Gt, Vb>& v0,
    const indexedVertex<Gt, Vb>& v1
)
{
    return v0.index_ == v1.type_ || v1.index_ == v0.type_;
}


template<class Gt, class Vb>
bool CGAL::boundaryTriangle
(
    const indexedVertex<Gt, Vb>& v0,
    const indexedVertex<Gt, Vb>& v1,
    const indexedVertex<Gt, Vb>& v2
)
{
    return (v0.pairPoint() && pointPair(v1, v2))
        || (v1.pairPoint() && pointPair(v2, v0))
        || (v2.pairPoint() && pointPair(v0, v1));
}


template<class Gt, class Vb>
bool CGAL::outsideTriangle
(
    const indexedVertex<Gt, Vb>& v0,
    const indexedVertex<Gt, Vb>& v1,
    const indexedVertex<Gt, Vb>& v2
)
{
    return (v0.farPoint() || v0.ppSlave())
        || (v1.farPoint() || v1.ppSlave())
        || (v2.farPoint() || v2.ppSlave());
}


// ************************************************************************* //
