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
    Foam::patchToPatch

Description
    Class to generate coupling geometry between two primitive patches

SourceFiles
    patchToPatch.C
    patchToPatchParallelOps.C

\*---------------------------------------------------------------------------*/

#ifndef patchToPatch_H
#define patchToPatch_H

#include "distributionMap.H"
#include "primitivePatch.H"
#include "primitiveOldTimePatch.H"
#include "remote.H"
#include "runTimeSelectionTables.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class patchToPatch Declaration
\*---------------------------------------------------------------------------*/

class patchToPatch
{
protected:

    // Private Data

        //- Flag to indicate that the two patches are co-directional and
        //  that the orientation of one should therefore be reversed
        const bool reverse_;

        //- Index of the processor holding all faces of the patchToPatch, or -1
        //  if spread across multiple processors
        label singleProcess_;

        //- For each source face, the coupled local target faces
        List<DynamicList<label>> srcLocalTgtFaces_;

        //- For each target face, the coupled local source faces
        List<DynamicList<label>> tgtLocalSrcFaces_;

        //- Map from source patch faces to target-local source patch faces
        autoPtr<distributionMap> srcMapPtr_;

        //- Map from target patch faces to source-local target patch faces
        autoPtr<distributionMap> tgtMapPtr_;

        //- When running in parallel, a map from local source face index to
        //  source processor and face index
        autoPtr<List<remote>> localSrcProcFacesPtr_;

        //- When running in parallel, a map from local target face index to
        //  target processor and face index
        autoPtr<List<remote>> localTgtProcFacesPtr_;


    // Private Member Functions

        // Searching

            //- Get the bound box for a source face
            virtual treeBoundBox srcBox
            (
                const face& srcFace,
                const pointField& srcPoints,
                const vectorField& srcPointNormals
            ) const = 0;

            //- Get the bound box for a source face
            treeBoundBox srcBox
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const label srcFacei
            ) const;

            //- Get the bound box for the source patch
            treeBoundBox srcBox
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0
            ) const;

            //- Get the bound box for a target face
            treeBoundBox tgtBox
            (
                const primitiveOldTimePatch& tgtPatch,
                const label tgtFacei
            ) const;

            //- Get the bound box for the target patch
            treeBoundBox tgtBox
            (
                const primitiveOldTimePatch& tgtPatch
            ) const;


        // Intersection

            //- Intersect two faces
            virtual bool intersectFaces
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const label srcFacei,
                const label tgtFacei
            ) = 0;

            //- Intersect two faces
            bool findOrIntersectFaces
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const label srcFacei,
                const label tgtFacei
            );

            //- Intersect a queue of source-target face pairs. Update completion
            //  lists and form a new queue of target-source face pairs.
            label intersectPatchQueue
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const bool isSrc,
                const DynamicList<labelPair>& queue,
                labelList& faceComplete,
                DynamicList<labelPair>& otherQueue,
                const labelList& otherFaceComplete,
                boolList& otherFaceQueued,
                boolList& otherFaceVisited
            );

            //- Intersect the patches
            void intersectPatches
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            );


        // Parallel functionality

            //- Determine which target faces need to be sent to the source.
            //  This is done before intersection. Bound boxes are used to
            //  estimate what faces will intersect.
            labelListList tgtPatchSendFaces
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            ) const;

            //- Distribute a patch given its distribution map
            static List<remote> distributePatch
            (
                const distributionMap& map,
                const primitiveOldTimePatch& patch,
                autoPtr<PrimitiveOldTimePatch<faceList, pointField>>&
                    localPatchPtr
            );


        // Hooks

            //- Initialisation
            virtual void initialise
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            );

            //- Trim the local target patch so that only parts that actually
            //  intersect the source remain
            virtual labelList trimLocalTgt
            (
                const primitiveOldTimePatch& localTgtPatch
            );

            //- Distribute the source patch so that any information needed by
            //  the target is locally available
            virtual void distributeSrc(const primitiveOldTimePatch& srcPatch);

            //- Send data that resulted from an intersection between the source
            //  patch and a distributed source-local-target patch back to the
            //  original target processes
            virtual void rDistributeTgt(const primitiveOldTimePatch& tgtPatch);

            //- Finalising
            virtual label finalise
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const transformer& tgtToSrc
            );


        // Interpolation

            //- For each source face, the coupled target weights
            virtual tmpNrc<List<DynamicList<scalar>>> srcWeights() const = 0;

            //- For each target face, the coupled source weights
            virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights() const = 0;

            //- Return a primitive with all components set to NaN
            template<class Type>
            static Type NaN();

            //- Interpolate with normalisation
            template<class Type>
            static tmp<Field<Type>> interpolate
            (
                const List<DynamicList<label>>& localOtherFaces,
                const List<DynamicList<scalar>>& weights,
                const autoPtr<distributionMap>& otherMapPtr,
                const Field<Type>& otherFld
            );

            //- Interpolate without normalisation
            template<class Type>
            static tmp<Field<Type>> interpolate
            (
                const List<DynamicList<label>>& localOtherFaces,
                const List<DynamicList<scalar>>& weights,
                const autoPtr<distributionMap>& otherMapPtr,
                const Field<Type>& otherFld,
                const Field<Type>& leftOverFld
            );


public:

    //- Runtime type information
    TypeName("patchToPatch");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        patchToPatch,
        bool,
        (const bool reverse),
        (reverse)
    );


    // Constructors

        //- Construct from components
        patchToPatch(const bool reverse);

        //- Disallow default bitwise copy construction
        patchToPatch(const patchToPatch&) = delete;


    //- Destructor
    virtual ~patchToPatch();


    // Selector

        //- Select from name
        static autoPtr<patchToPatch> New
        (
            const word& patchToPatchType,
            const bool reverse
        );


    // Member Functions

        // Access

            //- Flag to indicate that the two patches are co-directional and
            //  that the orientation of one should therefore be reversed
            inline bool reverse() const;

            //- Index of the processor holding all faces of the patchToPatch,
            //  or -1 if spread across multiple processors
            inline label singleProcess() const;

            //- Is this intersection on a single process?
            inline bool isSingleProcess() const;

            //- For each source face, the coupled target procs and faces
            inline List<List<remote>> srcTgtProcFaces() const;

            //- For each target face, the coupled source procs and faces
            inline List<List<remote>> tgtSrcProcFaces() const;


        // Interpolation

            //- Interpolate a source patch field to the target with no left
            //  over values specified. If the interpolation weight sum is less
            //  than one for a face then they will be normalised. If the
            //  interpolation weight sum is zero for a face then that face's
            //  value will be NaN.
            template<class Type>
            tmp<Field<Type>> srcToTgt(const Field<Type>& srcFld) const;

            //- Interpolate a source patch field to the target with left over
            //  values specified. If the interpolation weight sum is less than
            //  one for a face then the average will include the left over
            //  value multiplied by one minus the weight sum.
            template<class Type>
            tmp<Field<Type>> srcToTgt
            (
                const Field<Type>& srcFld,
                const Field<Type>& leftOverTgtFld
            ) const;

            //- Interpolate a target patch field to the source with no left
            //  over values specified. As the corresponding srcToTgt.
            template<class Type>
            tmp<Field<Type>> tgtToSrc(const Field<Type>& tgtFld) const;

            //- Interpolate a target patch field to the source with left
            //  over values specified. As the corresponding srcToTgt.
            template<class Type>
            tmp<Field<Type>> tgtToSrc
            (
                const Field<Type>& tgtFld,
                const Field<Type>& leftOverSrcFld
            ) const;


        // Manipulation

            //- Update addressing and weights for the given patches
            void update
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const transformer& tgtToSrc = NullObjectRef<transformer>()
            );

            //- Update addressing and weights for the given patches
            void update
            (
                const primitivePatch& srcPatch,
                const vectorField& srcPointNormals,
                const primitivePatch& tgtPatch,
                const transformer& tgtToSrc = NullObjectRef<transformer>()
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const patchToPatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchToPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchToPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
