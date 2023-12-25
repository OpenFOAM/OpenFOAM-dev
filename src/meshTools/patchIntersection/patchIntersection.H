/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::patchIntersection

Description
    Base class for patch intersections. Provides type name and debugging. See
    templated derivatives for actual functionality.

SourceFiles
    patchIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef patchIntersection_H
#define patchIntersection_H

#include "DynamicField.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class patchIntersection Declaration
\*---------------------------------------------------------------------------*/

class patchIntersection
{
protected:

    // Protected static data

        //- Flag to set whether the faces of the intersection are oriented the
        //  same way as the source patch (true) or the target patch (false).
        static const bool orientToSource_;


    // Protected Data

        // Points

            //- The intersection points
            DynamicField<point> points_;

            //- The source points' corresponding intersection points
            labelList srcPointPoints_;

            //- The target points' corresponding intersection points
            labelList tgtPointPoints_;

            //- The intersection points' corresponding source points, or -1
            //  if the point does not correspond to a source point
            DynamicList<label> pointSrcPoints_;

            //- The intersection points' corresponding target points, or -1
            //  if the point does not correspond to a target point
            DynamicList<label> pointTgtPoints_;

            //- The source edges' intersection points. Ordered along the edge.
            List<DynamicList<label>> srcEdgePoints_;

            //- The target edges' intersection points. Ordered along the edge.
            List<DynamicList<label>> tgtEdgePoints_;

            //- The intersection points' source edges, or -1 if the point
            //  is not on a source edge
            DynamicList<label> pointSrcEdges_;

            //- The intersection points' target edges, or -1 if the point
            //  is not on a target edge
            DynamicList<label> pointTgtEdges_;

            // !!! We don't store srcFacePoints and tgtFacePoints. These are
            // not needed at present. They also don't have any ordering
            // associated with them, so they don't need to be maintained by the
            // intersection process. They could be generated on demand as a
            // post-processing step from from pointSrcFaces and pointTgtFaces.

            //- The intersection points' source faces, or -1 if the point
            //  did not project into a source face
            DynamicList<label> pointSrcFaces_;

            //- The intersection points' target faces, or -1 if the point
            //  did not project into a target face
            DynamicList<label> pointTgtFaces_;


        // Faces

            //- The intersection faces
            DynamicList<face> faces_;

            //- The source faces' intersection faces
            List<DynamicList<label>> srcFaceFaces_;

            //- The target faces' intersection faces
            List<DynamicList<label>> tgtFaceFaces_;

            //- The intersection faces' corresponding source faces, or -1
            //  if the face does not correspond to source face
            DynamicList<label> faceSrcFaces_;

            //- The intersection faces' corresponding target faces, or -1
            //  if the face does not correspond to target face
            DynamicList<label> faceTgtFaces_;


    // Protected Member Functions

        //- Report properties of the intersection process
        void report(const word& writeSuffix=word::null);


public:

    // Runtime type information

        ClassName("patchIntersection");

        virtual word type() const = 0;


    // Constructors

        //- Construct given sizes
        patchIntersection
        (
            const label srcNPoints,
            const label tgtNPoints,
            const label srcNEdges,
            const label tgtNEdges,
            const label srcNFaces,
            const label tgtNFaces
        );


    //- Destructor
    virtual ~patchIntersection();


    // Member Functions

        // Points

            //- The intersection points
            inline const pointField& points() const;

            //- The source points' corresponding intersection points
            inline const labelList& srcPointPoints() const;

            //- The target points' corresponding intersection points
            inline const labelList& tgtPointPoints() const;

            //- The intersection points' corresponding source points, or -1
            //  if the point does not correspond to a source point
            inline const DynamicList<label>& pointSrcPoints() const;

            //- The intersection points' corresponding target points, or -1
            //  if the point does not correspond to a target point
            inline const DynamicList<label>& pointTgtPoints() const;

            //- The source edges' intersection points. Ordered along the edge.
            inline const List<DynamicList<label>>& srcEdgePoints() const;

            //- The target edges' intersection points. Ordered along the edge.
            inline const List<DynamicList<label>>& tgtEdgePoints() const;

            //- The intersection points' source edges, or -1 if the point
            //  is not on a source edge
            inline const DynamicList<label>& pointSrcEdges() const;

            //- The intersection points' target edges, or -1 if the point
            //  is not on a target edge
            inline const DynamicList<label>& pointTgtEdges() const;

            //- The intersection points' source faces, or -1 if the point
            //  did not project into a source face
            inline const DynamicList<label>& pointSrcFaces() const;

            //- The intersection points' target faces, or -1 if the point
            //  did not project into a target face
            inline const DynamicList<label>& pointTgtFaces() const;


        // Faces

            //- The intersection faces
            inline const faceList& faces() const;

            //- The source faces' intersection faces
            inline const List<DynamicList<label>>& srcFaceFaces() const;

            //- The target faces' intersection faces
            inline const List<DynamicList<label>>& tgtFaceFaces() const;

            //- The intersection faces' corresponding source faces, or -1
            //  if the face does not correspond to source face
            inline const DynamicList<label>& faceSrcFaces() const;

            //- The intersection faces' corresponding target faces, or -1
            //  if the face does not correspond to target face
            inline const DynamicList<label>& faceTgtFaces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchIntersectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
