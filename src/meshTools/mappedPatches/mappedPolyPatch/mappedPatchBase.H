/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::mappedPatchBase

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

    If constructed from dictionary:
    \verbatim
        // Region to sample (default is region0)
        sampleRegion region0;

        // What to sample:
        // - nearestCell         : sample cell containing point
        // - nearestOnlyCell     : nearest sample cell (even if not containing
        //                         point)
        // - nearestPatchFace    : nearest face on selected patch
        // - nearestPatchFaceAMI : nearest face on selected patch
                                   - patches need not conform
                                   - uses AMI interpolation
        // - nearestFace         : nearest boundary face on any patch
        sampleMode nearestCell;

        // If sampleMode is nearestPatchFace : patch to find faces of
        samplePatch movingWall;

        // If sampleMode is nearestPatchFace : specify patchgroup to find
        // samplePatch and sampleRegion (if not provided)
        coupleGroup baffleGroup;

        // How to supply offset (w.r.t. my patch face centres):
        // - uniform : single offset vector
        // - nonuniform : per-face offset vector
        // - normal : using supplied distance and face normal
        offsetMode uniform;

        // According to offsetMode (see above) supply one of
        // offset, offsets or distance
        offset  (1 0 0);
    \endverbatim

    Note: if offsetMode is \c normal it uses outwards pointing normals. So
    supply a negative distance if sampling inside the domain.

SourceFiles
    mappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchBase_H
#define mappedPatchBase_H

#include "pointField.H"
#include "AMIInterpolation.H"
#include "coupleGroupIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;
class distributionMap;

/*---------------------------------------------------------------------------*\
                       Class mappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class mappedPatchBase
{

public:

    // Type enumerations

        //- Mesh items to sample
        enum sampleMode
        {
            NEARESTCELL,         // nearest cell containing sample
            NEARESTPATCHFACE,    // nearest face on selected patch
            NEARESTPATCHFACEAMI, // nearest patch face + AMI interpolation
            NEARESTFACE          // nearest face
        };

        //- How to project face centres
        enum offsetMode
        {
            UNIFORM,            // single offset vector
            NONUNIFORM,         // per-face offset vector
            NORMAL              // use face normal + distance
        };

        static const NamedEnum<sampleMode, 4> sampleModeNames_;

        static const NamedEnum<offsetMode, 3> offsetModeNames_;


protected:

    // Protected data

        //- Patch to sample
        const polyPatch& patch_;

        //- Region to sample
        mutable word sampleRegion_;

        //- Same region
        mutable bool sameRegion_;

        //- What to sample
        const sampleMode mode_;

        //- Patch (if in sampleMode NEARESTPATCH*)
        mutable word samplePatch_;

        //- PatchGroup (if in sampleMode NEARESTPATCH*)
        const coupleGroupIdentifier coupleGroup_;

        //- How to obtain samples
        offsetMode offsetMode_;

        //- Offset vector (uniform)
        vector offset_;

        //- Offset vector (nonuniform)
        vectorField offsets_;

        //- Offset distance (normal)
        scalar distance_;


        // Derived information

            //- Distributor
            mutable autoPtr<distributionMap> mapPtr_;

            //- Map pre-addressing
            mutable labelList mapIndices_;


        // AMI interpolator (only for NEARESTPATCHFACEAMI)

            //- Pointer to AMI interpolator
            mutable autoPtr<AMIInterpolation> AMIPtr_;

            //- Flag to indicate that slave patch should be reversed for AMI
            const bool AMIReverse_;

            //- Pointer to projection surface employed by AMI interpolator
            mutable autoPtr<searchableSurface> surfPtr_;

            //- Dictionary storing projection surface description
            dictionary surfDict_;


    // Protected Member Functions

        //- Find cells/faces containing samples
        void findSamples
        (
            const globalIndex& patchGlobalIndex,
            labelList& sampleGlobalPatchFaces,
            labelList& sampleIndices
        ) const;

        //- Return size of mapped mesh/patch/boundary
        label sampleSize() const;

        //- Calculate mapping
        void calcMapping() const;

        //- Calculate AMI interpolator
        void calcAMI() const;

        //- Helper to read field or non-uniform list from dictionary
        static tmp<pointField> readListOrField
        (
            const word& keyword,
            const dictionary& dict,
            const label size
        );


public:

    //- Runtime type information
    TypeName("mappedPatchBase");


    // Constructors

        //- Construct from patch
        mappedPatchBase(const polyPatch&);

        //- Construct with offsetMode=non-uniform
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const vectorField& offsets
        );

        //- Construct from offsetMode=uniform
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const vector& offset
        );

        //- Construct from offsetMode=normal and distance
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar distance
        );

        //- Construct from dictionary
        mappedPatchBase(const polyPatch&, const dictionary&);

        //- Construct as copy, resetting patch
        mappedPatchBase(const polyPatch&, const mappedPatchBase&);

        //- Construct as copy, resetting patch, map original data
        mappedPatchBase
        (
            const polyPatch&,
            const mappedPatchBase&,
            const labelUList& mapAddressing
        );


    //- Destructor
    virtual ~mappedPatchBase();


    // Member Functions

        // Access

            //- What to sample
            inline const sampleMode& mode() const;

            //- Region to sample
            inline const word& sampleRegion() const;

            //- Patch (only if NEARESTPATCHFACE)
            inline const word& samplePatch() const;

            //- Cached sampleRegion != mesh.name()
            inline bool sameRegion() const;

            //- Return reference to the parallel distribution map
            inline const distributionMap& map() const;

            //- Return reference to the indices that have to be supplied to the
            //  parallel distribution map
            inline const labelList& mapIndices() const;

            //- Get the region mesh
            const polyMesh& sampleMesh() const;

            //- Get the patch on the region
            const polyPatch& samplePolyPatch() const;

            //- Get the sample points
            tmp<pointField> samplePoints() const;



        // Edit

            //- Clear out data on mesh change
            void clearOut();


        // Distribute

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            tmp<Field<Type>> distribute(const Field<Type>& fld) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            tmp<Field<Type>> distribute(const tmp<Field<Type>>& fld) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            tmp<Field<Type>> reverseDistribute(const Field<Type>& fld) const;

            //- Wrapper around map/interpolate data distribution
            template<class Type>
            tmp<Field<Type>> reverseDistribute
            (
                const tmp<Field<Type>>& fld
            ) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            tmp<Field<Type>> reverseDistribute
            (
                const Field<Type>& fld,
                const CombineOp& cop
            ) const;

            //- Wrapper around map/interpolate data distribution with operation
            template<class Type, class CombineOp>
            tmp<Field<Type>> reverseDistribute
            (
                const tmp<Field<Type>>& fld,
                const CombineOp& cop
            ) const;


        // I/O

            //- Write as a dictionary
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mappedPatchBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedPatchBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
