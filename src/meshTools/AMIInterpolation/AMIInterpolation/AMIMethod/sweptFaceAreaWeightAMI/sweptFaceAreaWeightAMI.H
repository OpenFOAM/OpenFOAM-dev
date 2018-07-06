/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::sweptFaceAreaWeightAMI

Description
    Swept face area weighted Arbitrary Mesh Interface (AMI) method

    This method uses the point normals of the source patch to sweep the source
    faces over the target patches. This creates a projection which fills space,
    and which therefore generates overlap areas which are consistent between
    neighbouring faces. The projection of a source edge is shown below:
    \verbatim
                /
               /
          n_1 ^
             /
            /
       p_1 o
           \\
            \\  ^ u
             \\  \
              \\  + - > v
               \\
                o - - - > - - -
              p_0      n_0
    \endverbatim

    The surface is, in general, not flat. Any deviation between the two end
    normals generates expansion, contraction, and twist in the surface. The
    surface connects with surfaces emanating from connected edges along the end
    normals. This is what makes the projection fill space and generate
    consistent overlaps between neighbouring faces.

    The projected surface is parameterised by the local coordinates, \f$u\f$
    and \f$v\f$, and a position on this plane is calculated from \f$u\f$ and
    \f$v\f$ as follows:
    \f[
        x(u, v) = p_0 + (p_1 - p_0) u + [ q_0 + (q_1 - q_0) u ] v
    \f]

    To calculate an intersection with a line between points \f$l_0\f$ to
    \f$l_1\f$, we define a local coordinate, \f$w\f$, along the line, and
    subtract it's equation from that of the surface:
    \f[
        0 = (p_0 - l_0) - (l_1 - l_0) w + (p_1 - p_0) u +
            [ q_0 + (q_1 - q_0) u ] v
    \f]

    This is a system of three equations in three unknowns. It is non-linear,
    courtesy of the \f$u v\f$ term at the end. It can be reduced to a single
    quadratic in any of the three variables. We choose to solve for \f$u\f$ by
    taking the dot product of the above equation with the following vector:
    \f[
        (l_1 - l_0) \times [ q_0 + (q_1 - q_0) u ]
    \f]

    The sign of the intersection (i.e., whether the line crosses from below the
    surface to above or vice versa) can be determined by taking the dot product
    of the line vector with the surface normal at the intersection. The surface
    normal is as follows:
    \f[
        n(u, v) = (p_1 - p_0) \times [q_0 + (q_1 - q_0) u] + (q_1 \times q_0) v
    \f]

SourceFiles
    sweptFaceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef sweptFaceAreaWeightAMI_H
#define sweptFaceAreaWeightAMI_H

#include "faceAreaWeightAMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sweptFaceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class sweptFaceAreaWeightAMI
:
    public faceAreaWeightAMI<SourcePatch, TargetPatch>
{
private:

    // Private static data

        //- Minimum length of a cut as a ratio of the overall projected edge
        //  length
        static const scalar minCutRatio_;

        //- Maximum allowable dot product between a source point normal and a
        //  target triangle
        static const scalar maxDot_;


    // Private classes

        //- A fixed list of tris which simulates a dynamic list by incrementing
        //  a counter whenever its append method is called. This is used as an
        //  optimisation so that the tri cutting does not allocate memory.
        template<unsigned Size>
        class cutTriList
        :
            public FixedList<FixedList<point, 3>, Size>
        {
        private:

            //- The number of stored elements
            label n_;


        public:

            //- Construct null
            cutTriList()
            :
                n_(0)
            {}

            //- Clear the array
            void clear()
            {
                n_ = 0;
            }

            //- Get the current size
            label size() const
            {
                return n_;
            }

            //- Add a new tet to the end of the array
            void append(const FixedList<point, 3>& t)
            {
                this->operator[](n_) = t;
                ++ n_;
            }
        };


    // Private member functions

        // Debugging

            //- Write a VTK file of cut triangles
            template<unsigned Size>
            void writeCutTrisVTK
            (
                const cutTriList<Size>& tris,
                const word& name
            ) const;

            //- Write an OBJ file of a face
            void writeFaceOBJ
            (
                const face& f,
                const pointField& ps,
                const string& name
            ) const;

            //- Write an OBJ file of the source projection
            void writeProjectionOBJ
            (
                const label srcN,
                const FixedList<point, 4>& srcTri,
                const FixedList<point, 4>& srcPrj
            ) const;


        //- Convert the source tris and normals to a projection. Most of the
        //  time this does nothing, but if some of the normals point in the
        //  reverse direction the projection will be reduced to span only the
        //  region in which the projection points forward through the target
        //  plane. Returns the number of edges in the projection (0, 3 or 4).
        label getSourceProjection
        (
            FixedList<point, 4>& srcTri,
            FixedList<point, 4>& srcNrm,
            const FixedList<point, 3>& tgtTri
        ) const;

        //- Get the cutting plane, for an edge of the source projection.
        plane getCutPlane
        (
            const point& p0,
            const point& p1,
            const vector& n0,
            const vector& n1,
            const FixedList<point, 3>& tgtTri
        ) const;

        //- The minimum weight below which connections are discarded
        virtual scalar minWeight() const;

        //- The maximum edge angle that the walk will cross
        virtual scalar maxWalkAngle() const;



protected:

    // Protected Member Functions

        // Evaluation

            //- Area of intersection between source and target faces
            virtual scalar interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("sweptFaceAreaWeightAMI");


    // Constructors

        //- Use parent constructors
        using faceAreaWeightAMI<SourcePatch, TargetPatch>::faceAreaWeightAMI;


    //- Destructor
    virtual ~sweptFaceAreaWeightAMI();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sweptFaceAreaWeightAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
