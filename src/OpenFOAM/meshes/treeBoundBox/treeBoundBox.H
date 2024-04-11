/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::treeBoundBox

Description
    Standard boundBox + extra functionality for use in octree.

    Numbering of corner points is according to octant numbering.

    On the back plane (z=0):

    \verbatim
        Y
        ^
        |
        +--------+
        |2      3|
        |        |
        |        |
        |        |
        |0      1|
        +--------+->X
    \endverbatim

    For the front plane add 4 to the point labels.


SourceFiles
    treeBoundBoxI.H
    treeBoundBox.C
    treeBoundBoxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef treeBoundBox_H
#define treeBoundBox_H

#include "boundBox.H"
#include "direction.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class randomGenerator;


// Forward declaration of friend functions and operators

class treeBoundBox;

bool operator==(const treeBoundBox&, const treeBoundBox&);
bool operator!=(const treeBoundBox&, const treeBoundBox&);

Istream& operator>>(Istream& is, treeBoundBox&);
Ostream& operator<<(Ostream& os, const treeBoundBox&);


/*---------------------------------------------------------------------------*\
                        Class treeBoundBox Declaration
\*---------------------------------------------------------------------------*/

class treeBoundBox
:
    public boundBox
{
        //- To initialise edges.
        static edgeList calcEdges(const label[12][2]);

        //- To initialise faceNormals.
        static FixedList<vector, 6> calcFaceNormals();

public:

    // Static Data Members

        //- The great value used for greatBox and invertedBox
        static const scalar great;

        //- As per boundBox::greatBox, but with great instead of vGreat
        static const treeBoundBox greatBox;

        //- As per boundBox::invertedBox, but with great instead of vGreat
        static const treeBoundBox invertedBox;

        //- Bits used for octant/point coding.
        //  Every octant/corner point is the combination of three faces.
        struct octantBit
        {
            enum
            {
                rightHalf = 0x1 << 0,
                topHalf   = 0x1 << 1,
                frontHalf = 0x1 << 2
            };
        };

        //- Face codes
        struct faceId
        {
            enum
            {
                left   = 0,
                right  = 1,
                bottom = 2,
                top    = 3,
                back   = 4,
                front  = 5
            };
        };

        //- Bits used for face coding
        struct faceBit
        {
            enum
            {
                noFace    = 0,
                left      = 0x1 << faceId::left,    //1
                right     = 0x1 << faceId::right,   //2
                bottom    = 0x1 << faceId::bottom,  //4
                top       = 0x1 << faceId::top,     //8
                back      = 0x1 << faceId::back,    //16
                front     = 0x1 << faceId::front,   //32
            };
        };

        //- Edges codes.
        //  e01 = edge between 0 and 1.
        struct edgeId
        {
            enum
            {
                e01 = 0,
                e13 = 1,
                e23 = 2,
                e02 = 3,

                e45 = 4,
                e57 = 5,
                e67 = 6,
                e46 = 7,

                e04 = 8,
                e15 = 9,
                e37 = 10,
                e26 = 11
            };
        };

        //- Face to point addressing
        static const faceList faces;

        //- Edge to point addressing
        static const edgeList edges;

        //- Per face the unit normal
        static const FixedList<vector, 6> faceNormals;


    // Constructors

        //- Construct null setting points to zero
        inline treeBoundBox();

        //- Construct from a boundBox
        explicit inline treeBoundBox(const boundBox& bb);

        //- Construct from components
        inline treeBoundBox(const point& min, const point& max);

        //- Construct as the bounding box of the given pointField.
        //  Local processor domain only (no reduce as in boundBox)
        explicit treeBoundBox(const UList<point>&);

        //- Construct as subset of points
        //  Local processor domain only (no reduce as in boundBox)
        treeBoundBox(const UList<point>&, const labelUList& indices);

        //- Construct as subset of points
        //  The indices could be from edge/triFace etc.
        //  Local processor domain only (no reduce as in boundBox)
        template<unsigned Size>
        treeBoundBox
        (
            const UList<point>&,
            const FixedList<label, Size>& indices
        );

        //- Construct from dictionary
        inline treeBoundBox(const dictionary&);

        //- Construct from Istream
        inline treeBoundBox(Istream&);


    // Member Functions

        // Access

            //- Typical dimension length,height,width
            inline scalar typDim() const;

            //- Vertex coordinates. In octant coding.
            tmp<pointField> points() const;


        // Check

            //- Corner point given octant
            inline point corner(const direction) const;

            //- Sub box given by octant number. Midpoint calculated.
            treeBoundBox subBbox(const direction) const;

            //- Sub box given by octant number. Midpoint provided.
            treeBoundBox subBbox(const point& mid, const direction) const;

            //- Returns octant number given point and the calculated midpoint.
            inline direction subOctant
            (
                const point& pt
            ) const;

            //- Returns octant number given point and midpoint.
            static inline direction subOctant
            (
                const point& mid,
                const point& pt
            );

            //- Returns octant number given point and the calculated midpoint.
            //  onEdge set if the point is on edge of subOctant
            inline direction subOctant
            (
                const point& pt,
                bool& onEdge
            ) const;

            //- Returns octant number given point and midpoint.
            //  onEdge set if the point is on edge of subOctant
            static inline direction subOctant
            (
                const point& mid,
                const point& pt,
                bool& onEdge
            );

            //- Returns octant number given intersection and midpoint.
            //  onEdge set if the point is on edge of subOctant
            //  If onEdge, the direction vector determines which octant to use
            //  (acc. to which octant the point would be if it were moved
            //  along dir)
            static inline direction subOctant
            (
                const point& mid,
                const vector& dir,
                const point& pt,
                bool& onEdge
            );

            //- Calculates optimal order to look for nearest to point.
            //  First will be the octant containing the point,
            //  second the octant with boundary nearest to the point etc.
            inline void searchOrder
            (
                 const point& pt,
                 FixedList<direction, 8>& octantOrder
            ) const;

            //- Overlaps other bounding box?
            using boundBox::overlaps;

            //- Intersects segment; set point to intersection position and face,
            //  return true if intersection found.
            //  (pt argument used during calculation even if not intersecting).
            //  Calculates intersections from outside supplied vector
            //  (overallStart, overallVec). This is so when
            //  e.g. tracking through lots of consecutive boxes
            // (typical octree) we're not accumulating truncation errors. Set
            // to start, (end-start) if not used.
            bool intersects
            (
                const point& overallStart,
                const vector& overallVec,
                const point& start,
                const point& end,
                point& pt,
                direction& ptBits
            ) const;

            //- Like above but does not return faces point is on
            bool intersects
            (
                const point& start,
                const point& end,
                point& pt
            ) const;

            //- Contains point or other bounding box?
            using boundBox::contains;

            //- Contains point (inside or on edge) and moving in direction
            //  dir would cause it to go inside.
            bool contains(const vector& dir, const point&) const;

            //- Code position of point on bounding box faces
            direction faceBits(const point&) const;

            //- Position of point relative to bounding box
            direction posBits(const point&) const;

            //- Calculate nearest and furthest (to point) vertex coords of
            //  bounding box
            void calcExtremities
            (
                const point& pt,
                point& nearest,
                point& furthest
            ) const;

            //- Returns distance point to furthest away corner.
            scalar maxDist(const point&) const;

            //- Compare distance to point with other bounding box
            //  return:
            //  -1 : all vertices of my bounding box are nearer than any of
            //       other
            //  +1 : all vertices of my bounding box are further away than
            //       any of other
            //   0 : none of the above.
            label distanceCmp(const point&, const treeBoundBox& other) const;

            //- Return asymetrically extended bounding box, with guaranteed
            //  minimum width of s*mag(span) in any direction
            inline treeBoundBox extend(const scalar s) const;

    // Write

        void writeOBJ(const fileName& fName) const;


    // Friend Operators

        friend bool operator==(const treeBoundBox&, const treeBoundBox&);
        friend bool operator!=(const treeBoundBox&, const treeBoundBox&);


    // IOstream operator

        friend Istream& operator>>(Istream& is, treeBoundBox&);
        friend Ostream& operator<<(Ostream& os, const treeBoundBox&);
};


//- Data associated with treeBoundBox type are contiguous
template<>
inline bool contiguous<treeBoundBox>() {return contiguous<boundBox>();}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "treeBoundBoxI.H"

#ifdef NoRepository
    #include "treeBoundBoxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
