/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "Pstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex()
:
    Vb(),
    type_(vtUnassigned),
    index_(-1),
    processor_(Foam::Pstream::myProcNo()),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(const Point& p)
:
    Vb(p),
    type_(vtUnassigned),
    index_(-1),
    processor_(Foam::Pstream::myProcNo()),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex
(
    const Point& p,
    vertexType type
)
:
    Vb(p),
    type_(type),
    index_(-1),
    processor_(Foam::Pstream::myProcNo()),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex
(
    const Foam::point& p,
    vertexType type
)
:
    Vb(Point(p.x(), p.y(), p.z())),
    type_(type),
    index_(-1),
    processor_(Foam::Pstream::myProcNo()),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex
(
    const Point& p,
    Foam::label index,
    vertexType type,
    int processor
)
:
    Vb(p),
    type_(type),
    index_(index),
    processor_(processor),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex
(
    const Foam::point& p,
    Foam::label index,
    vertexType type,
    int processor
)
:
    Vb(Point(p.x(), p.y(), p.z())),
    type_(type),
    index_(index),
    processor_(processor),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(const Point& p, Cell_handle f)
:
    Vb(f, p),
    type_(vtUnassigned),
    index_(-1),
    processor_(Foam::Pstream::myProcNo()),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(Cell_handle f)
:
    Vb(f),
    type_(vtUnassigned),
    index_(-1),
    processor_(Foam::Pstream::myProcNo()),
    alignment_(Foam::triad::unset),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Gt, class Vb>
inline Foam::label& CGAL::indexedVertex<Gt, Vb>::index()
{
    return index_;
}


template<class Gt, class Vb>
inline Foam::label CGAL::indexedVertex<Gt, Vb>::index() const
{
    return index_;
}


template<class Gt, class Vb>
inline Foam::indexedVertexEnum::vertexType&
CGAL::indexedVertex<Gt, Vb>::type()
{
    return type_;
}


template<class Gt, class Vb>
inline Foam::indexedVertexEnum::vertexType
CGAL::indexedVertex<Gt, Vb>::type() const
{
    return type_;
}


template<class Gt, class Vb>
inline Foam::tensor& CGAL::indexedVertex<Gt, Vb>::alignment()
{
    return alignment_;
}


template<class Gt, class Vb>
inline const Foam::tensor& CGAL::indexedVertex<Gt, Vb>::alignment() const
{
    return alignment_;
}


template<class Gt, class Vb>
inline Foam::scalar& CGAL::indexedVertex<Gt, Vb>::targetCellSize()
{
    return targetCellSize_;
}


template<class Gt, class Vb>
inline Foam::scalar CGAL::indexedVertex<Gt, Vb>::targetCellSize() const
{
    return targetCellSize_;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::farPoint() const
{
    return type_ == vtFar;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalPoint() const
{
    return type_ == vtInternal || type_ == vtInternalNearBoundary;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::referred() const
{
    // Can't be zero as the first few points are far points which won't be
    // referred
    // return index_ < 0;

    // processor_ will be take the value of the processor that this vertex is
    // from, so it cannot be on this processor.
    return processor_ != Foam::Pstream::myProcNo();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::real() const
{
    return (internalPoint() || boundaryPoint()) && !referred();
}


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::procIndex() const
{
    return processor_;
}


template<class Gt, class Vb>
inline int& CGAL::indexedVertex<Gt, Vb>::procIndex()
{
    return processor_;
}


template<class Gt, class Vb>
inline void CGAL::indexedVertex<Gt, Vb>::setInternal()
{
    type_ = vtInternal;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::nearBoundary() const
{
    return type_ == vtInternalNearBoundary;
}


template<class Gt, class Vb>
inline void CGAL::indexedVertex<Gt, Vb>::setNearBoundary()
{
    type_ = vtInternalNearBoundary;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::boundaryPoint() const
{
    return type_ >= vtInternalSurface && !farPoint();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalOrBoundaryPoint() const
{
    return internalPoint() || internalBoundaryPoint();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::nearOrOnBoundary() const
{
    return boundaryPoint() || nearBoundary();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalBoundaryPoint() const
{
    return type_ >= vtInternalSurface && type_ <= vtInternalFeaturePoint;
}

template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalBaffleSurfacePoint() const
{
    return (type_ == vtInternalSurfaceBaffle);
}

template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalBaffleEdgePoint() const
{
    return (type_ == vtInternalFeatureEdgeBaffle);
}

template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::externalBoundaryPoint() const
{
    return type_ >= vtExternalSurface && type_ <= vtExternalFeaturePoint;
}

template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::externalBaffleSurfacePoint() const
{
    return (type_ == vtExternalSurfaceBaffle);
}

template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::externalBaffleEdgePoint() const
{
    return (type_ == vtExternalFeatureEdgeBaffle);
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::constrained() const
{
    return type_ == vtConstrained;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::featurePoint() const
{
    return type_ == vtInternalFeaturePoint || type_ == vtExternalFeaturePoint;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::featureEdgePoint() const
{
    return type_ == vtInternalFeatureEdge || type_ == vtExternalFeatureEdge;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::surfacePoint() const
{
    return type_ == vtInternalSurface || type_ == vtExternalSurface;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::fixed() const
{
    return vertexFixed_;
}


template<class Gt, class Vb>
inline bool& CGAL::indexedVertex<Gt, Vb>::fixed()
{
    return vertexFixed_;
}


// ************************************************************************* //
