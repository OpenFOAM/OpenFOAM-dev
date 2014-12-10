/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "PatchTools.H"
#include "ListListOps.H"
#include "mergePoints.H"
#include "face.H"
#include "triFace.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PatchTools::gatherAndMerge
(
    const scalar mergeDist,
    const PrimitivePatch<Face, FaceList, PointField, PointType>& p,
    Field<PointType>& mergedPoints,
    List<Face>& mergedFaces,
    labelList& pointMergeMap
)
{
    // Collect points from all processors
    labelList pointSizes;
    {
        List<Field<PointType> > gatheredPoints(Pstream::nProcs());
        gatheredPoints[Pstream::myProcNo()] = p.points();

        Pstream::gatherList(gatheredPoints);

        if (Pstream::master())
        {
            pointSizes = ListListOps::subSizes
            (
                gatheredPoints,
                accessOp<Field<PointType> >()
            );

            mergedPoints = ListListOps::combine<Field<PointType> >
            (
                gatheredPoints,
                accessOp<Field<PointType> >()
            );
        }
    }

    // Collect faces from all processors and renumber using sizes of
    // gathered points
    {
        List<List<Face> > gatheredFaces(Pstream::nProcs());
        gatheredFaces[Pstream::myProcNo()] = p;
        Pstream::gatherList(gatheredFaces);

        if (Pstream::master())
        {
            mergedFaces = static_cast<const List<Face>&>
            (
                ListListOps::combineOffset<List<Face> >
                (
                    gatheredFaces,
                    pointSizes,
                    accessOp<List<Face> >(),
                    offsetOp<Face>()
                )
            );
        }
    }

    if (Pstream::master())
    {
        Field<PointType> newPoints;
        labelList oldToNew;

        bool hasMerged = mergePoints
        (
            mergedPoints,
            mergeDist,
            false,                  // verbosity
            oldToNew,
            newPoints
        );

        if (hasMerged)
        {
            // Store point mapping
            pointMergeMap.transfer(oldToNew);

            // Copy points
            mergedPoints.transfer(newPoints);

            // Relabel faces
            List<Face>& faces = mergedFaces;

            forAll(faces, faceI)
            {
                inplaceRenumber(pointMergeMap, faces[faceI]);
            }
        }
    }
}


// ************************************************************************* //
