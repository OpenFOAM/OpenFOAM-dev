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
public:

    // Public Structures

        //- Structure to conveniently store processor and face indices
        struct procFace
        {
            //- The processor index
            label proci;

            //- The face index
            label facei;

            //- Equality comparison
            friend bool operator==(const procFace& a, const procFace& b)
            {
                return a.proci == b.proci && a.facei == b.facei;
            }

            //- Inequality comparison
            friend bool operator!=(const procFace& a, const procFace& b)
            {
                return !(a == b);
            }

            //- Output stream operator
            friend Ostream& operator<<(Ostream& os, const procFace& p)
            {
                return os << p.proci << token::SPACE << p.facei;
            }

            //- Input stream operator
            friend Istream& operator>>(Istream& is, procFace& p)
            {
                return is >> p.proci >> p.facei;
            }
        };


protected:

    // Private Data

        //- Flag to indicate that the two patches are co-directional and
        //  that the orientation of one should therefore be reversed
        const bool reverse_;

        //- Index of the processor holding all faces of the patchToPatch, or -1
        //  if spread across multiple processors
        label singleProcess_;

        //- When running in parallel, a map from local source face index to
        //  source processor and face index
        autoPtr<List<procFace>> localSrcProcFacesPtr_;

        //- When running in parallel, a map from local target face index to
        //  target processor and face index
        autoPtr<List<procFace>> localTgtProcFacesPtr_;

        //- For each source face, the coupled local target faces
        List<DynamicList<label>> srcLocalTgtFaces_;

        //- For each target face, the coupled local source faces
        List<DynamicList<label>> tgtLocalSrcFaces_;


    // Private Member Functions

        // Indexing

            //- Transfer list-list b into list-list a
            template<class SubListA, class SubListB>
            static inline void transferListList
            (
                List<SubListA>& a,
                List<SubListB>& b
            );

            //- Reverse distribute a list-list given the map
            template<class Type>
            static inline void rDistributeListList
            (
                const label size,
                const distributionMap& map,
                List<List<Type>>& data
            );

            //- Reverse distribute a dynamic list-list given the map
            template<class Type>
            static inline void rDistributeListList
            (
                const label size,
                const distributionMap& map,
                List<DynamicList<Type>>& data
            );

            //- Map local faces to proc faces
            static List<List<procFace>> localFacesToProcFaces
            (
                const List<DynamicList<label>>& localFaces,
                const List<procFace>& map = NullObjectRef<List<procFace>>()
            );

            //- Map proc faces to local faces
            static List<DynamicList<label>> procFacesToLocalFaces
            (
                const List<List<procFace>>& procFaces,
                const HashTable<label, procFace, Hash<procFace>>& map
            );


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

            //- Determine whether or not the intersection of the given patches
            //  is on a single process. Set singleProcess_ to that process if
            //  so, and to -1 if not.
            void calcSingleProcess
            (
                const primitiveOldTimePatch& srcPatch,
                const primitiveOldTimePatch& tgtPatch
            );

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

            //- Determine which source faces need to be sent to the target.
            //  This is done after intersection, using the addressing generated
            //  by the intersection.
            labelListList srcPatchSendFaces() const;

            //- Create a distribution map from a list-list of faces to be sent
            //  (i.e., the result of calling one of the above methods).
            distributionMap patchDistributionMap
            (
                labelListList&& sendFaces
            ) const;

            //- Distribute a patch given its distribution map. Just generate
            //  local-proc-face addressing; not the distributed patch itself.
            void distributePatch
            (
                const distributionMap& map,
                List<procFace>& localProcFaces
            ) const;

            //- Distribute a patch given its distribution map
            PrimitiveOldTimePatch<faceList, pointField> distributePatch
            (
                const distributionMap& map,
                const primitiveOldTimePatch& patch,
                List<procFace>& localProcFaces
            ) const;


        // Hooks

            //- Initialisation
            virtual void initialise
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch
            );

            //- Distribute the target patch so that enough is locally available
            //  for its intersection with the source patch can be computed.
            //  Happens before intersection. Bound boxes are used to determine
            //  what is needed where. Return the resulting local patch and the
            //  map from which it was calculated.
            virtual
            tmpNrc<PrimitiveOldTimePatch<faceList, pointField>>
            distributeTgt
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                distributionMap& tgtMap
            );

            //- Distribute the source patch so that everything the target
            //  intersects is locally available. Happens after intersection.
            //  The intersection addressing is used to determine what is needed
            //  where. Return the resulting local patch and the map from which
            //  it was calculated.
            virtual
            tmpNrc<PrimitiveOldTimePatch<faceList, pointField>>
            distributeSrc
            (
                const primitiveOldTimePatch& srcPatch,
                distributionMap& srcMap
            );

            //- Send data that resulted from an intersection between the source
            //  patch and a distributed source-local-target patch back to the
            //  original target processes.
            virtual void rDistributeTgt
            (
                const primitiveOldTimePatch& tgtPatch,
                const distributionMap& tgtMap
            );

            //- Finalising
            virtual label finalise
            (
                const primitiveOldTimePatch& srcPatch,
                const vectorField& srcPointNormals,
                const vectorField& srcPointNormals0,
                const primitiveOldTimePatch& tgtPatch,
                const transformer& tgtToSrc
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
            inline List<List<procFace>> srcTgtProcFaces() const;

            //- For each target face, the coupled source procs and faces
            inline List<List<procFace>> tgtSrcProcFaces() const;

            //- For each source face, the coupled target weights
            virtual tmpNrc<List<DynamicList<scalar>>> srcWeights
            (
                const primitivePatch& srcPatch
            ) const = 0;

            //- For each target face, the coupled source weights
            virtual tmpNrc<List<DynamicList<scalar>>> tgtWeights
            (
                const primitivePatch& tgtPatch
            ) const = 0;


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

#endif

// ************************************************************************* //
