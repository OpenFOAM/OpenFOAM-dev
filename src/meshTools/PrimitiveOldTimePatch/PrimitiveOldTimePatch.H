/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::PrimitiveOldTimePatch

Description

SourceFiles
    PrimitiveOldTimePatch.C

\*---------------------------------------------------------------------------*/

#ifndef PrimitiveOldTimePatch_H
#define PrimitiveOldTimePatch_H

#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Trait to determine what face list should be used for the old-time patch
template<class Container>
struct UFaceList;
template<class Face>
struct UFaceList<UList<Face>> { typedef UList<Face> type; };
template<class Face>
struct UFaceList<List<Face>> : public UFaceList<UList<Face>> {};
template<class Face>
struct UFaceList<SubList<Face>> : public UFaceList<UList<Face>> {};
template<class Face>
struct UFaceList<UIndirectList<Face>> { typedef UIndirectList<Face> type; };
template<class Face>
struct UFaceList<IndirectList<Face>> : public UFaceList<UIndirectList<Face>> {};

/*---------------------------------------------------------------------------*\
                     Class PrimitiveOldTimePatch Declaration
\*---------------------------------------------------------------------------*/

template<class FaceList, class PointField>
class PrimitiveOldTimePatch
:
    public PrimitivePatch<FaceList, PointField>
{
public:

    // Public Typedefs

        using typename PrimitivePatch<FaceList, PointField>::FaceListType;

        using typename PrimitivePatch<FaceList, PointField>::FaceType;

        using typename PrimitivePatch<FaceList, PointField>::PointFieldType;

        using typename PrimitivePatch<FaceList, PointField>::PointType;


private:

    // Private Typedefs

        typedef
            PrimitivePatch
            <
                typename UFaceList<FaceListType>::type,
                const Field<PointType>&
            >
            patch0Type;


    // Private Data

        //- Stored old-time points
        autoPtr<Field<PointType>> points0Ptr_;

        //- Reference to global list of old-time points
        const Field<PointType>& points0_;

        //- Engine for calculating old-time geometry. Note: Methods that
        //  generate topology should not be called here. The base patch should
        //  have all the necessary topology available. The calculation of some
        //  geometric quantities requires topological data to be available. In
        //  these cases, special steps need to be taken here to make sure that
        //  the old-time patch object does not generate duplicate topology.
        autoPtr<patch0Type> patch0Ptr_;

        //- Points local to patch
        mutable Field<PointType>* localPoints0Ptr_;


public:

    // Constructors

        //- Construct from components
        PrimitiveOldTimePatch
        (
            const FaceList& faces,
            const Field<PointType>& points,
            const Field<PointType>& points0
        );

        //- Construct from patch and old-time points
        PrimitiveOldTimePatch
        (
            const PrimitivePatch<FaceList, PointField>& patch,
            const Field<PointType>& points0
        );

        //- Construct from components without old-time points
        PrimitiveOldTimePatch
        (
            const FaceList& faces,
            const Field<PointType>& points
        );

        //- Construct from patch without old-time points
        PrimitiveOldTimePatch
        (
            const PrimitivePatch<FaceList, PointField>& patch
        );

        //- Copy constructor
        PrimitiveOldTimePatch(const PrimitiveOldTimePatch& patch);

        //- Move constructor
        PrimitiveOldTimePatch(PrimitiveOldTimePatch&& patch);

        //- Construct and return a clone
        autoPtr<PrimitiveOldTimePatch<FaceList, PointField>> clone() const
        {
            return autoPtr<PrimitiveOldTimePatch<FaceList, PointField>>
            (
                new PrimitiveOldTimePatch<FaceList, PointField>(*this)
            );
        }


    //- Destructor
    virtual ~PrimitiveOldTimePatch();


    // Member Functions

        // Access

            //- Return whether or not old-time geometry is available
            bool has0() const
            {
                return patch0Ptr_.valid();
            }

            //- Return reference to old-time global points
            const Field<PointType>& points0() const
            {
                return points0_;
            }


        // Geometry

            //- Return pointField of old-time points in patch
            const Field<PointType>& localPoints0() const;

            //- Return old-time face centres for patch
            const Field<PointType>& faceCentres0() const;

            //- Return old-time face areas for patch
            const Field<PointType>& faceAreas0() const;

            //- Return old-time face normals for patch
            const Field<PointType>& faceNormals0() const;

            //- Return old-time point normals for patch
            const Field<PointType>& pointNormals0() const;


        // Edit

            //- ...
            void clearOut();

            //- ...
            void clearGeom();

            //- Correct patch after moving points
            virtual void movePoints0(const Field<PointType>&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PrimitiveOldTimePatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PrimitiveOldTimePatch.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
