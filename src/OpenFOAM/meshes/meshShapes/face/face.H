/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::face

Description
    A face is a list of labels corresponding to mesh vertices.

See also
    Foam::triFace

SourceFiles
    faceI.H
    face.C
    faceIntersection.C
    faceContactSphere.C
    faceAreaInContact.C
    faceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef face_H
#define face_H

#include "pointField.H"
#include "labelList.H"
#include "edgeList.H"
#include "vectorField.H"
#include "faceListFwd.H"
#include "intersection.H"
#include "pointHit.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class face;
class triFace;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicList;

inline bool operator==(const face& a, const face& b);
inline bool operator!=(const face& a, const face& b);
inline Istream& operator>>(Istream&, face&);

/*---------------------------------------------------------------------------*\
                            Class face Declaration
\*---------------------------------------------------------------------------*/

class face
:
    public labelList
{
public:

    //- Return types for classify
    enum proxType
    {
        NONE,
        POINT,  // Close to point
        EDGE    // Close to edge
    };

    // Static Data Members

        static const char* const typeName;


    // Static Functions

        //- Return centre point given face points
        template<class PointField>
        static vector centre(const PointField& ps);

        //- Return vector area given face points
        template<class PointField>
        static vector area(const PointField& ps);


    // Constructors

        //- Construct null
        inline face();

        //- Construct given size
        explicit inline face(label);

        //- Construct from list of labels
        explicit inline face(const labelUList&);

        //- Construct from list of labels
        explicit inline face(const labelList&);

        //- Construct by transferring the parameter contents
        explicit inline face(labelList&&);

        //- Copy construct from triFace
        face(const triFace&);

        //- Construct from Istream
        inline face(Istream&);


    // Member Functions

        //- Collapse face by removing duplicate point labels
        //  return the collapsed size
        label collapse();

        //- Flip the face in-place.
        //  The starting points of the original and flipped face are identical.
        void flip();

        //- Return the points corresponding to this face
        inline pointField points(const pointField&) const;

        //- Centre point of face
        point centre(const pointField&) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average(const pointField&, const Field<Type>&) const;

        //- Return vector area
        vector area(const pointField&) const;

        //- Return scalar magnitude
        inline scalar mag(const pointField&) const;

        //- Return unit normal
        vector normal(const pointField&) const;

        //- Return face with reverse direction
        //  The starting points of the original and reverse face are identical.
        face reverseFace() const;

        //- Navigation through face vertices

            //- Which vertex on face (face index given a global index)
            //  returns -1 if not found
            label which(const label globalIndex) const;

            //- Next vertex on face
            inline label nextLabel(const label i) const;

            //- Previous vertex on face
            inline label prevLabel(const label i) const;


        //- Return the volume swept out by the face when its points move
        scalar sweptVol
        (
            const pointField& oldPoints,
            const pointField& newPoints
        ) const;

        //- Return the inertia tensor, with optional reference
        //  point and density specification
        tensor inertia
        (
            const pointField&,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return potential intersection with face with a ray starting
        //  at p, direction n (does not need to be normalised)
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. Face-centre is calculated from point average.
        //  For a hit, the distance is signed.  Positive number
        //  represents the point in front of triangle
        //  In case of miss the point is the nearest point on the face
        //  and the distance is the distance between the intersection point
        //  and the original point.
        //  The half-ray or full-ray intersection and the contact
        //  sphere adjustment of the projection vector is set by the
        //  intersection parameters
        pointHit ray
        (
            const point& p,
            const vector& n,
            const pointField&,
            const intersection::algorithm alg =
                intersection::algorithm::fullRay,
            const intersection::direction dir =
                intersection::direction::vector
        ) const;

        //- Fast intersection with a ray.
        //  Does face-centre decomposition and returns triangle intersection
        //  point closest to p. See triangle::intersection for details.
        pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const pointField&,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        pointHit nearestPoint
        (
            const point& p,
            const pointField&
        ) const;

        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2)
        //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
        //    Note: edges are counted from starting vertex so
        //    e.g. edge n is from f[n] to f[0], where the face has n + 1
        //    points
        pointHit nearestPointClassify
        (
            const point& p,
            const pointField&,
            label& nearType,
            label& nearLabel
        ) const;

        //- Return contact sphere diameter
        scalar contactSphereDiameter
        (
            const point& p,
            const vector& n,
            const pointField&
        ) const;

        //- Return area in contact, given the displacement in vertices
        scalar areaInContact
        (
            const pointField&,
            const scalarField& v
        ) const;

        //- Return number of edges
        inline label nEdges() const;

        //- Return edges in face point ordering,
        //  i.e. edges()[0] is edge between [0] and [1]
        edgeList edges() const;

        //- Return n-th face edge
        inline edge faceEdge(const label n) const;

        //- Return the edge direction on the face
        //  Returns:
        //  -  0: edge not found on the face
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        int edgeDirection(const edge&) const;

        //- Size of the face's triangulation
        inline label nTriangles() const;

        //- Compare faces
        //   0: different
        //  +1: identical
        //  -1: same face, but different orientation
        static int compare(const face&, const face&);

        //- Return true if the faces have the same vertices
        static bool sameVertices(const face&, const face&);


    // Member Operators

        //- Move assignment labelList
        inline void operator=(labelList&&);


    // Friend Operators

        friend bool operator==(const face& a, const face& b);
        friend bool operator!=(const face& a, const face& b);


    // Istream Operator

        friend Istream& operator>>(Istream&, face&);
};


//- Hash specialisation to offset faces in ListListOps::combineOffset
template<>
class offsetOp<face>
{

public:

    inline face operator()
    (
        const face& x,
        const label offset
    ) const
    {
        face result(x.size());

        forAll(x, xI)
        {
            result[xI] = x[xI] + offset;
        }
        return result;
    }
};


// Global functions

//- Find the longest edge on a face. Face point labels index into pts.
label longestEdge(const face& f, const pointField& pts);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceI.H"

#ifdef NoRepository
    #include "faceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
