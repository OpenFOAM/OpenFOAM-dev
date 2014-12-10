/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
PrimitivePatch
(
    const FaceList<Face>& faces,
    const Field<PointType>& points
)
:
    FaceList<Face>(faces),
    points_(points),
    edgesPtr_(NULL),
    nInternalEdges_(-1),
    boundaryPointsPtr_(NULL),
    faceFacesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    pointEdgesPtr_(NULL),
    pointFacesPtr_(NULL),
    localFacesPtr_(NULL),
    meshPointsPtr_(NULL),
    meshPointMapPtr_(NULL),
    edgeLoopsPtr_(NULL),
    localPointsPtr_(NULL),
    localPointOrderPtr_(NULL),
    faceCentresPtr_(NULL),
    faceNormalsPtr_(NULL),
    pointNormalsPtr_(NULL)
{}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
PrimitivePatch
(
    const Xfer<FaceList<Face> >& faces,
    const Xfer<List<PointType> >& points
)
:
    FaceList<Face>(faces),
    points_(points),
    edgesPtr_(NULL),
    nInternalEdges_(-1),
    boundaryPointsPtr_(NULL),
    faceFacesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    pointEdgesPtr_(NULL),
    pointFacesPtr_(NULL),
    localFacesPtr_(NULL),
    meshPointsPtr_(NULL),
    meshPointMapPtr_(NULL),
    edgeLoopsPtr_(NULL),
    localPointsPtr_(NULL),
    localPointOrderPtr_(NULL),
    faceCentresPtr_(NULL),
    faceNormalsPtr_(NULL),
    pointNormalsPtr_(NULL)
{}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
PrimitivePatch
(
    FaceList<Face>& faces,
    Field<PointType>& points,
    const bool reUse
)
:
    FaceList<Face>(faces, reUse),
    points_(points, reUse),
    edgesPtr_(NULL),
    nInternalEdges_(-1),
    boundaryPointsPtr_(NULL),
    faceFacesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    pointEdgesPtr_(NULL),
    pointFacesPtr_(NULL),
    localFacesPtr_(NULL),
    meshPointsPtr_(NULL),
    meshPointMapPtr_(NULL),
    edgeLoopsPtr_(NULL),
    localPointsPtr_(NULL),
    localPointOrderPtr_(NULL),
    faceCentresPtr_(NULL),
    faceNormalsPtr_(NULL),
    pointNormalsPtr_(NULL)
{}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
PrimitivePatch
(
    const PrimitivePatch<Face, FaceList, PointField, PointType>& pp
)
:
    PrimitivePatchName(),
    FaceList<Face>(pp),
    points_(pp.points_),
    edgesPtr_(NULL),
    nInternalEdges_(-1),
    boundaryPointsPtr_(NULL),
    faceFacesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    pointEdgesPtr_(NULL),
    pointFacesPtr_(NULL),
    localFacesPtr_(NULL),
    meshPointsPtr_(NULL),
    meshPointMapPtr_(NULL),
    edgeLoopsPtr_(NULL),
    localPointsPtr_(NULL),
    localPointOrderPtr_(NULL),
    faceCentresPtr_(NULL),
    faceNormalsPtr_(NULL),
    pointNormalsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
~PrimitivePatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
movePoints
(
    const Field<PointType>&
)
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "movePoints() : "
            << "recalculating PrimitivePatch geometry following mesh motion"
            << endl;
    }

    clearGeom();
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::edgeList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
edges() const
{
    if (!edgesPtr_)
    {
        calcAddressing();
    }

    return *edgesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::label
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
nInternalEdges() const
{
    if (!edgesPtr_)
    {
        calcAddressing();
    }

    return nInternalEdges_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
boundaryPoints() const
{
    if (!boundaryPointsPtr_)
    {
        calcBdryPoints();
    }

    return *boundaryPointsPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelListList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
faceFaces() const
{
    if (!faceFacesPtr_)
    {
        calcAddressing();
    }

    return *faceFacesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelListList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        calcAddressing();
    }

    return *edgeFacesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelListList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        calcAddressing();
    }

    return *faceEdgesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelListList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }

    return *pointEdgesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelListList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
pointFaces() const
{
    if (!pointFacesPtr_)
    {
        calcPointFaces();
    }

    return *pointFacesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::List<Face>&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
localFaces() const
{
    if (!localFacesPtr_)
    {
        calcMeshData();
    }

    return *localFacesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
meshPoints() const
{
    if (!meshPointsPtr_)
    {
        calcMeshData();
    }

    return *meshPointsPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::Map<Foam::label>&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
meshPointMap() const
{
    if (!meshPointMapPtr_)
    {
        calcMeshPointMap();
    }

    return *meshPointMapPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::Field<PointType>&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
localPoints() const
{
    if (!localPointsPtr_)
    {
        calcLocalPoints();
    }

    return *localPointsPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
localPointOrder() const
{
    if (!localPointOrderPtr_)
    {
        calcLocalPointOrder();
    }

    return *localPointOrderPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::label
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
whichPoint
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


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::Field<PointType>&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
faceCentres() const
{
    if (!faceCentresPtr_)
    {
        calcFaceCentres();
    }

    return *faceCentresPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::Field<PointType>&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
faceNormals() const
{
    if (!faceNormalsPtr_)
    {
        calcFaceNormals();
    }

    return *faceNormalsPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::Field<PointType>&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        calcPointNormals();
    }

    return *pointNormalsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
operator=
(
    const PrimitivePatch<Face, FaceList, PointField, PointType>& pp
)
{
    clearOut();

    FaceList<Face>::operator=(pp);
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
