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

    Example:
    \verbatim
        // What to sample
        // - nearestCell         : cell containing the sampling point
        // - nearestPatchFace    : nearest face on selected patch
        // - nearestPatchFaceAMI : AMI interpolation from selected patch
        // - nearestFace         : nearest boundary face on any patch
        sampleMode nearestCell;

        // Region to sample from
        sampleRegion region0;

        // If sampleMode is nearestPatchFace* : The patch to sample from
        samplePatch movingWall;

        // If sampleMode is nearestPatchFace* : Couple group to specify the
        //                                      sample region and patch
        //                                      (alternative to specifying
        //                                      sampleRegion and samplePatch
        //                                      directly)
        coupleGroup baffleGroup;

        // How to offset the patch face centres to the sampling locations
        // - none                : no offset
        // - normal              : distance along the patch face normals
        // - direction           : specified offset vector
        offsetMode direction;

        // If offsetMode is normal : The normal distance to offset
        distance 1;

        // If offsetMode is direction : The offset vector
        offset (1 0 0);
    \endverbatim

    Note that patch normals point outward, so if offsetMode is \c normal then a
    negative distance will be required in order to sample value from inside the
    domain.

SourceFiles
    mappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPatchBase_H
#define mappedPatchBase_H

#include "pointField.H"
#include "AMIInterpolation.H"
#include "patchToPatch.H"
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
            PATCHTOPATCH,        // ...
            NEARESTFACE          // nearest face
        };

        //- How to project face centres
        enum offsetMode
        {
            NONE,       // do not offset
            NORMAL,     // offset a distance in the face-normal direction
            DIRECTION   // offset with a specified vector
        };

        static const NamedEnum<sampleMode, 5> sampleModeNames_;

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

        //- How to offset the patch face centres to the sampling locations
        offsetMode offsetMode_;

        //- Offset distance
        scalar distance_;

        //- Offset vector
        vector offset_;


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


        //- Patch-to-patch intersection engine (only for PATCHTOPATCH)

            //- Is the patch-to-patch intersection engine up to date?
            mutable bool patchToPatchIsValid_;

            //- Patch-to-patch intersection engine
            mutable autoPtr<patchToPatch> patchToPatchPtr_;


    // Protected Member Functions

        //- Read the offset mode from a dictionary
        offsetMode readOffsetMode(const dictionary& dict) const;

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

        //- Calculate the patch-to-patch intersection engine
        void calcPatchToPatch() const;

        //- Helper to read field or non-uniform list from dictionary
        static tmp<pointField> readListOrField
        (
            const word& keyword,
            const dictionary& dict,
            const label size
        );

        //- Return whether or not the sample patch (if any) is of mapped type
        inline bool sampleIsMappedPatch() const;

        //- Get the mapped sample patch
        inline const mappedPatchBase& sampleMappedPatch() const;


public:

    //- Runtime type information
    TypeName("mappedPatchBase");


    // Constructors

        //- Construct from patch
        mappedPatchBase(const polyPatch&);

        //- Construct with no offset
        mappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch
        );

        //- Construct from dictionary
        mappedPatchBase(const polyPatch&, const dictionary&);

        //- Construct as copy, resetting patch
        mappedPatchBase(const polyPatch&, const mappedPatchBase&);


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
            virtual tmp<pointField> samplePoints() const;



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
