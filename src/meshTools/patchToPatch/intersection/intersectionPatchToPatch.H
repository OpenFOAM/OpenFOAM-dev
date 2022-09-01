/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::patchToPatches::intersection

Description
    Class to generate patchToPatch coupling geometry. A full geometric
    intersection is done between a face and those opposite, and coupling
    geometry is calculated accordingly.

SourceFiles
    intersection.C

\*---------------------------------------------------------------------------*/

#ifndef intersectionPatchToPatch_H
#define intersectionPatchToPatch_H

#include "patchToPatch.H"
#include "polygonTriangulate.H"
#include "triFaceList.H"
#include "triIntersectLocation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{

/*---------------------------------------------------------------------------*\
                          Class intersection Declaration
\*---------------------------------------------------------------------------*/

class intersection
:
    public patchToPatch
{
public:

    // Public Structures

        //- Structure to store the geometry associated with part of a patch
        //  face
        struct part
        {
            //- The area of this part
            vector area;

            //- The centre of this part
            point centre;

            //- Default construct
            part()
            :
                area(Zero),
                centre(NaN<point>())
            {}

            //- Default construct
            part(const zero&)
            :
                part()
            {}

            //- Construct from an area and a centre
            part(const vector& a, const point& c)
            :
                area(a),
                centre(c)
            {}

            //- Construct from a polygon
            template<class PointField>
            part(const PointField& ps)
            :
                area(face::area(ps)),
                centre(face::centre(ps))
            {}

            //- Negate this part
            part operator-() const
            {
                return part(- area, centre);
            }

            //- Add another part to this one
            void operator+=(const part& p)
            {
                const scalar magArea = mag(area);
                const scalar magPArea = mag(p.area);

                area = area + p.area;

                centre =
                    magArea == 0 ? p.centre
                  : magPArea == 0 ? centre
                  : (magArea*centre + magPArea*p.centre)/(magArea + magPArea);
            }

            //- Subtract another part from this one
            void operator-=(const part& p)
            {
                this->operator+=(-p);
            }

            //- Equality comparison
            friend bool operator==(const part& a, const part& b)
            {
                return a.area == b.area && a.centre == b.centre;
            }

            //- Inequality comparison
            friend bool operator!=(const part& a, const part& b)
            {
                return !(a == b);
            }

            //- Output stream operator
            friend Ostream& operator<<(Ostream& os, const part& p)
            {
                return os << p.area << p.centre;
            }

            //- Input stream operator
            friend Istream& operator>>(Istream& is, part& p)
            {
                return is >> p.area >> p.centre;
            }
        };

        //- Structure to store the geometry associated with the coupling
        //  between parts of two patch faces
        struct couple
        :
            public part
        {
            //- The neighbour part
            part nbr;

            //- Default construct
            couple()
            :
                part(),
                nbr()
            {}

            //- Construct from a coupled pair of parts
            couple(const part& p, const part& nbrP)
            :
                part(p),
                nbr(nbrP)
            {}

            //- Equality comparison
            friend bool operator==(const couple& a, const couple& b)
            {
                return
                    static_cast<const part&>(a) == static_cast<const part&>(b)
                 && a.nbr == b.nbr;
            }

            //- Inequality comparison
            friend bool operator!=(const couple& a, const couple& b)
            {
                return !(a == b);
            }

            //- Output stream operator
            friend Ostream& operator<<(Ostream& os, const couple& c)
            {
                return os << static_cast<const part&>(c) << c.nbr;
            }

            //- Input stream operator
            friend Istream& operator>>(Istream& is, couple& c)
            {
                return is >> static_cast<part&>(c) >> c.nbr;
            }
        };


private:

    // Private Member Data

        // Geometry

            //- The coupling geometry for for each source face
            List<DynamicList<couple>> srcCouples_;

            //- The proportion of the source faces that are coupled
            List<scalar> srcCoverage_;

            //- The non-coupled geometry associated with each source edge
            List<part> srcEdgeParts_;

            //- The non-coupled geometry associated with mismatch in each
            //  source face's couplings
            List<part> srcErrorParts_;

            //- The coupling geometry for for each target face
            List<DynamicList<couple>> tgtCouples_;

            //- The proportion of the target faces that are coupled
            List<scalar> tgtCoverage_;


        //- Triangulation engine
        mutable polygonTriangulate triEngine_;


        // Workspace

            //- Source face triangulation points
            mutable List<triFaceList> srcTriPoints_;

            //- Source face triangulation edges
            mutable List<List<FixedList<label, 3>>> srcTriFaceEdges_;

            //- Target face triangulation points
            mutable List<triFaceList> tgtTriPoints_;

            //- Target face triangulation edges
            mutable List<List<FixedList<label, 3>>> tgtTriFaceEdges_;

            //- Source intersection points
            mutable DynamicList<point> ictSrcPoints_;

            //- Source intersection point normals
            mutable DynamicList<vector> ictSrcPointNormals_;

            //- Target intersection points
            mutable DynamicList<point> ictTgtPoints_;

            //- Intersection locations
            mutable DynamicList<triIntersect::location> ictPointLocations_;

            //- Source face edge parts
            DynamicList<part> srcFaceEdgePart_;

            //- Target face edge parts
            DynamicList<part> tgtFaceEdgePart_;

            //- Source face edge parts
            List<List<part>> srcFaceEdgeParts_;


        // Debugging

            //- Intersection points
            mutable DynamicList<point> debugPoints_;

            //- Intersection faces
            mutable DynamicList<labelList> debugFaces_;

            //- The source face associated with each intersection face
            mutable DynamicList<label> debugFaceSrcFaces_;

            //- The target face associated with each intersection face
            mutable DynamicList<label> debugFaceTgtFaces_;

            //- The side of the intersection each face is on; 1 is the source
            //  side, -1 is the target side, and 0 is a face that spans the
            //  intersection volume; e.g., an edge or error face.
            mutable DynamicList<label> debugFaceSides_;


    // Private Static Member Functions

        //- Return the values at the points of a tri face
        template<class Type>
        static FixedList<Type, 3> triPointValues
        (
            const triFace& t,
            const UList<Type>& values
        );


    // Private Member Functions

        //- Get the bound box for a source face
        virtual treeBoundBox srcBox
        (
            const face& srcFace,
            const pointField& srcPoints,
            const vectorField& srcPointNormals
        ) const;

        //- Intersect two faces
        virtual bool intersectFaces
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const label srcFacei,
            const label tgtFacei
        );

        //- Initialise the workspace
        virtual void initialise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch
        );

        //- Send data that resulted from an intersection between the source
        //  patch and a distributed source-local-target patch back to the
        //  original target processes.
        virtual void rDistributeTgt(const primitiveOldTimePatch& tgtPatch);

        //- Finalise the intersection
        virtual label finalise
        (
            const primitiveOldTimePatch& srcPatch,
            const vectorField& srcPointNormals,
            const vectorField& srcPointNormals0,
            const primitiveOldTimePatch& tgtPatch,
            const transformer& tgtToSrc
        );

        //- For each source face, the coupled target weights
        virtual tmpNrc<List<DynamicList<scalar>>> srcWeights() const;

        //- For each target face, the coupled source weights
        virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights() const;


public:

    //- Runtime type information
    TypeName("intersection");


    // Static Data Members

        //- Extra debugging for intersections between specific faces. Named
        //  "intersectionSrcFace" and "intersectionTgtFace" respectively.
        static int debugSrcFacei, debugTgtFacei;


    // Static Member Functions

        //- Get the bound box for a source face
        static treeBoundBox srcBoxStatic
        (
            const face& srcFace,
            const pointField& srcPoints,
            const vectorField& srcPointNormals
        );


    // Constructors

        //- Construct from components
        intersection(const bool reverse);


    //- Destructor
    ~intersection();


    // Member Functions

        //- For each source face, the source and target areas for each
        //  target coupling
        inline const List<DynamicList<couple>>& srcCouples() const;

        //- For each source edge, the non-coupled geometry associated
        //  with its projection
        inline const List<part>& srcEdgeParts() const;

        //- For each source face, the area associated with mismatch
        //  across the coupling
        inline const List<part>& srcErrorParts() const;

        //- For each target face, the target and source areas for each
        //  source coupling
        inline const List<DynamicList<couple>>& tgtCouples() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatches
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "intersectionPatchToPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
