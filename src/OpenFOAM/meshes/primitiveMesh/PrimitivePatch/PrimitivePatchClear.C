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

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "demandDrivenData.H"


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
clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    deleteDemandDrivenData(localPointsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearTopology()
{
    if (debug)
    {
        InfoInFunction << "Clearing patch addressing" << endl;
    }

    // group created and destroyed together
    if (edgesPtr_ && faceFacesPtr_ && edgeFacesPtr_ && faceEdgesPtr_)
    {
        delete edgesPtr_;
        edgesPtr_ = nullptr;

        delete faceFacesPtr_;
        faceFacesPtr_ = nullptr;

        delete edgeFacesPtr_;
        edgeFacesPtr_ = nullptr;

        delete faceEdgesPtr_;
        faceEdgesPtr_ = nullptr;
    }

    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(pointEdgesPtr_);
    deleteDemandDrivenData(pointFacesPtr_);
    deleteDemandDrivenData(edgeLoopsPtr_);
    deleteDemandDrivenData(localPointOrderPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearPatchMeshAddr()
{
    if (debug)
    {
        InfoInFunction << "Clearing patch-mesh addressing" << endl;
    }

    deleteDemandDrivenData(meshPointsPtr_);
    deleteDemandDrivenData(meshPointMapPtr_);
    deleteDemandDrivenData(localFacesPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
clearOut()
{
    clearGeom();
    clearTopology();
    clearPatchMeshAddr();
}


// ************************************************************************* //
