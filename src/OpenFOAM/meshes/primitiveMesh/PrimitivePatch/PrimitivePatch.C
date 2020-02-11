/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "Map.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    const FaceList& faces,
    const Field<PointType>& points
)
:
    FaceList(faces),
    points_(points),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    FaceList&& faces,
    Field<PointType>&& points
)
:
    FaceList(move(faces)),
    points_(move(points)),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    FaceList&& faces,
    List<PointType>&& points
)
:
    FaceList(move(faces)),
    points_(move(points)),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    FaceList& faces,
    Field<PointType>& points,
    const bool reuse
)
:
    FaceList(faces, reuse),
    points_(points, reuse),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    const PrimitivePatch<FaceList, PointField>& pp
)
:
    PrimitivePatchName(),
    FaceList(pp),
    points_(pp.points_),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    PrimitivePatch<FaceList, PointField>&& pp
)
:
    PrimitivePatchName(),
    FaceList(move(pp)),
    points_(move(pp.points_)),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::PrimitivePatch
(
    Istream& is,
    const Field<PointType>& points
)
:
    FaceList(is),
    points_(points),
    edgesPtr_(nullptr),
    nInternalEdges_(-1),
    boundaryPointsPtr_(nullptr),
    faceFacesPtr_(nullptr),
    edgeFacesPtr_(nullptr),
    faceEdgesPtr_(nullptr),
    pointEdgesPtr_(nullptr),
    pointFacesPtr_(nullptr),
    localFacesPtr_(nullptr),
    meshPointsPtr_(nullptr),
    meshPointMapPtr_(nullptr),
    edgeLoopsPtr_(nullptr),
    localPointsPtr_(nullptr),
    localPointOrderPtr_(nullptr),
    faceCentresPtr_(nullptr),
    faceAreasPtr_(nullptr),
    faceNormalsPtr_(nullptr),
    pointNormalsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::PrimitivePatch<FaceList, PointField>::~PrimitivePatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::movePoints
(
    const Field<PointType>&
)
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
            << "movePoints() : "
            << "recalculating PrimitivePatch geometry following mesh motion"
            << endl;
    }

    clearGeom();
}


template<class FaceList, class PointField>
const Foam::edgeList& Foam::PrimitivePatch<FaceList, PointField>::edges() const
{
    if (!edgesPtr_)
    {
        calcAddressing();
    }

    return *edgesPtr_;
}


template<class FaceList, class PointField>
Foam::label Foam::PrimitivePatch<FaceList, PointField>::nInternalEdges() const
{
    if (!edgesPtr_)
    {
        calcAddressing();
    }

    return nInternalEdges_;
}


template<class FaceList, class PointField>
const Foam::labelList&
Foam::PrimitivePatch<FaceList, PointField>::boundaryPoints() const
{
    if (!boundaryPointsPtr_)
    {
        calcBdryPoints();
    }

    return *boundaryPointsPtr_;
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::
faceFaces() const
{
    if (!faceFacesPtr_)
    {
        calcAddressing();
    }

    return *faceFacesPtr_;
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        calcAddressing();
    }

    return *edgeFacesPtr_;
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        calcAddressing();
    }

    return *faceEdgesPtr_;
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }

    return *pointEdgesPtr_;
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::pointFaces() const
{
    if (!pointFacesPtr_)
    {
        calcPointFaces();
    }

    return *pointFacesPtr_;
}


template<class FaceList, class PointField>
const Foam::List
<
    typename Foam::PrimitivePatch<FaceList, PointField>::FaceType
>&
Foam::PrimitivePatch<FaceList, PointField>::localFaces() const
{
    if (!localFacesPtr_)
    {
        calcMeshData();
    }

    return *localFacesPtr_;
}


template<class FaceList, class PointField>
const Foam::labelList&
Foam::PrimitivePatch<FaceList, PointField>::meshPoints() const
{
    if (!meshPointsPtr_)
    {
        calcMeshData();
    }

    return *meshPointsPtr_;
}


template<class FaceList, class PointField>
const Foam::Map<Foam::label>&
Foam::PrimitivePatch<FaceList, PointField>::meshPointMap() const
{
    if (!meshPointMapPtr_)
    {
        calcMeshPointMap();
    }

    return *meshPointMapPtr_;
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitivePatch<FaceList, PointField>::PointType
>&
Foam::PrimitivePatch<FaceList, PointField>::localPoints() const
{
    if (!localPointsPtr_)
    {
        calcLocalPoints();
    }

    return *localPointsPtr_;
}


template<class FaceList, class PointField>
const Foam::labelList&
Foam::PrimitivePatch<FaceList, PointField>::localPointOrder() const
{
    if (!localPointOrderPtr_)
    {
        calcLocalPointOrder();
    }

    return *localPointOrderPtr_;
}


template<class FaceList, class PointField>
Foam::label Foam::PrimitivePatch<FaceList, PointField>::whichPoint
(
    const label gp
) const
{
    Map<label>::const_iterator fnd = meshPointMap().find(gp);

    if (fnd != meshPointMap().end())
    {
        return fnd();
    }
    else
    {
        // Not found
        return -1;
    }
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitivePatch<FaceList, PointField>::PointType
>&
Foam::PrimitivePatch<FaceList, PointField>::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        calcFaceCentres();
    }

    return *faceCentresPtr_;
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitivePatch<FaceList, PointField>::PointType
>&
Foam::PrimitivePatch<FaceList, PointField>::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        calcFaceAreas();
    }

    return *faceAreasPtr_;
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitivePatch<FaceList, PointField>::PointType
>&
Foam::PrimitivePatch<FaceList, PointField>::faceNormals() const
{
    if (!faceNormalsPtr_)
    {
        calcFaceNormals();
    }

    return *faceNormalsPtr_;
}


template<class FaceList, class PointField>
const Foam::Field
<
    typename Foam::PrimitivePatch<FaceList, PointField>::PointType
>&
Foam::PrimitivePatch<FaceList, PointField>::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        calcPointNormals();
    }

    return *pointNormalsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::operator=
(
    const PrimitivePatch<FaceList, PointField>& pp
)
{
    clearOut();

    FaceList::shallowCopy(pp);
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::operator=
(
    PrimitivePatch<FaceList, PointField>&& pp
)
{
    clearOut();

    FaceList::operator=(move(pp));

    // This is only valid if PointField is not a reference
    // points_ = move(pp.points_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PrimitivePatchAddressing.C"
#include "PrimitivePatchEdgeLoops.C"
#include "PrimitivePatchClear.C"
#include "PrimitivePatchBdryPoints.C"
#include "PrimitivePatchLocalPointOrder.C"
#include "PrimitivePatchMeshData.C"
#include "PrimitivePatchMeshEdges.C"
#include "PrimitivePatchPointAddressing.C"
#include "PrimitivePatchProjectPoints.C"
#include "PrimitivePatchCheck.C"

// ************************************************************************* //
