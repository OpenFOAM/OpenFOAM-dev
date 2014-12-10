/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "Map.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcMeshData() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcMeshData() : "
               "calculating mesh data in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate meshPoints
    // if they have already been calculated.
    if (meshPointsPtr_ || localFacesPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcMeshData()"
        )   << "meshPointsPtr_ or localFacesPtr_already allocated"
            << abort(FatalError);
    }

    // Create a map for marking points.  Estimated size is 4 times the
    // number of faces in the patch
    Map<label> markedPoints(4*this->size());


    // Important:
    // ~~~~~~~~~~
    // In <= 1.5 the meshPoints would be in increasing order but this gives
    // problems in processor point synchronisation where we have to find out
    // how the opposite side would have allocated points.

    ////- 1.5 code:
    //// if the point is used, set the mark to 1
    //forAll(*this, facei)
    //{
    //    const Face& curPoints = this->operator[](facei);
    //
    //    forAll(curPoints, pointi)
    //    {
    //        markedPoints.insert(curPoints[pointi], -1);
    //    }
    //}
    //
    //// Create the storage and store the meshPoints.  Mesh points are
    //// the ones marked by the usage loop above
    //meshPointsPtr_ = new labelList(markedPoints.toc());
    //labelList& pointPatch = *meshPointsPtr_;
    //
    //// Sort the list to preserve compatibility with the old ordering
    //sort(pointPatch);
    //
    //// For every point in map give it its label in mesh points
    //forAll(pointPatch, pointi)
    //{
    //    markedPoints.find(pointPatch[pointi])() = pointi;
    //}

    //- Unsorted version:
    DynamicList<label> meshPoints(2*this->size());
    forAll(*this, facei)
    {
        const Face& curPoints = this->operator[](facei);

        forAll(curPoints, pointi)
        {
            if (markedPoints.insert(curPoints[pointi], meshPoints.size()))
            {
                meshPoints.append(curPoints[pointi]);
            }
        }
    }
    // Transfer to straight list (reuses storage)
    meshPointsPtr_ = new labelList(meshPoints, true);


    // Create local faces. Note that we start off from copy of original face
    // list (even though vertices are overwritten below). This is done so
    // additional data gets copied (e.g. region number of labelledTri)
    localFacesPtr_ = new List<Face>(*this);
    List<Face>& lf = *localFacesPtr_;

    forAll(*this, facei)
    {
        const Face& curFace = this->operator[](facei);
        lf[facei].setSize(curFace.size());

        forAll(curFace, labelI)
        {
            lf[facei][labelI] = markedPoints.find(curFace[labelI])();
        }
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcMeshData() : "
               "finished calculating mesh data in PrimitivePatch"
            << endl;
    }
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
calcMeshPointMap() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcMeshPointMap() : "
               "calculating mesh point map in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate meshPoints
    // if they have already been calculated.
    if (meshPointMapPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcMeshPointMap()"
        )   << "meshPointMapPtr_ already allocated"
            << abort(FatalError);
    }

    const labelList& mp = meshPoints();

    meshPointMapPtr_ = new Map<label>(2*mp.size());
    Map<label>& mpMap = *meshPointMapPtr_;

    forAll(mp, i)
    {
        mpMap.insert(mp[i], i);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcMeshPointMap() : "
               "finished calculating mesh point map in PrimitivePatch"
            << endl;
    }
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
calcLocalPoints() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcLocalPoints() : "
               "calculating localPoints in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate localPoints
    // if they have already been calculated.
    if (localPointsPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcLocalPoints()"
        )   << "localPointsPtr_already allocated"
            << abort(FatalError);
    }

    const labelList& meshPts = meshPoints();

    localPointsPtr_ = new Field<PointType>(meshPts.size());

    Field<PointType>& locPts = *localPointsPtr_;

    forAll(meshPts, pointi)
    {
        locPts[pointi] = points_[meshPts[pointi]];
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            << "calcLocalPoints() : "
            << "finished calculating localPoints in PrimitivePatch"
            << endl;
    }
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
calcPointNormals() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcPointNormals() : "
               "calculating pointNormals in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate pointNormals
    // if they have already been calculated.
    if (pointNormalsPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcPointNormals()"
        )   << "pointNormalsPtr_already allocated"
            << abort(FatalError);
    }

    const Field<PointType>& faceUnitNormals = faceNormals();

    const labelListList& pf = pointFaces();

    pointNormalsPtr_ = new Field<PointType>
    (
        meshPoints().size(),
        PointType::zero
    );

    Field<PointType>& n = *pointNormalsPtr_;

    forAll(pf, pointi)
    {
        PointType& curNormal = n[pointi];

        const labelList& curFaces = pf[pointi];

        forAll(curFaces, facei)
        {
            curNormal += faceUnitNormals[curFaces[facei]];
        }

        curNormal /= mag(curNormal) + VSMALL;
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcPointNormals() : "
               "finished calculating pointNormals in PrimitivePatch"
            << endl;
    }
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
calcFaceCentres() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcFaceCentres() : "
               "calculating faceCentres in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (faceCentresPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcFaceCentres()"
        )   << "faceCentresPtr_already allocated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new Field<PointType>(this->size());

    Field<PointType>& c = *faceCentresPtr_;

    forAll(c, facei)
    {
        c[facei] = this->operator[](facei).centre(points_);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcFaceCentres() : "
               "finished calculating faceCentres in PrimitivePatch"
            << endl;
    }
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
calcFaceNormals() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcFaceNormals() : "
               "calculating faceNormals in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceNormals
    // if they have already been calculated.
    if (faceNormalsPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatch<Face, FaceList, PointField, PointType>::"
            "calcFaceNormals()"
        )   << "faceNormalsPtr_already allocated"
            << abort(FatalError);
    }

    faceNormalsPtr_ = new Field<PointType>(this->size());

    Field<PointType>& n = *faceNormalsPtr_;

    forAll(n, facei)
    {
        n[facei] = this->operator[](facei).normal(points_);
        n[facei] /= mag(n[facei]) + VSMALL;
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
               "calcFaceNormals() : "
               "finished calculating faceNormals in PrimitivePatch"
            << endl;
    }
}


// ************************************************************************* //
