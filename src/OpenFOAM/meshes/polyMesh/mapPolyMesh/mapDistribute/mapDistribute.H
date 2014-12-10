/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::mapDistribute

Description
    Class containing processor-to-processor mapping information.

    We store mapping from the bits-to-send to the complete starting list
    (subXXXMap) and from the received bits to their location in the new
    list (constructXXXMap).

Note:
    Schedule is a list of processor pairs (one send, one receive. One of
    them will be myself) which forms a scheduled (i.e. non-buffered) exchange.
    See distribute on how to use it.
    Note2: number of items sent on one processor have to equal the number
    of items received on the other processor.

    To aid constructing these maps there are the constructors from global
    numbering, either with or without transforms.

    - without transforms:
    Constructors using compact numbering: layout is
    - all my own elements first (whether used or not)
    - followed by used-only remote elements sorted by remote processor.
    So e.g 4 procs and on proc 1 the compact
    table will first have all globalIndex.localSize() elements from proc1
    followed by used-only elements of proc0, proc2, proc3.
    The constructed mapDistribute sends the local elements from and
    receives the remote elements into their compact position.
    compactMap[procI] is the position of elements from procI in the compact
    map. compactMap[myProcNo()] is empty since trivial addressing.

    It rewrites the input global indices into indices into the constructed
    data.


    - with transforms:
    This requires the precalculated set of possible transforms
    (globalIndexAndTransform). These are given as permutations (+, -, or none)
    of up to 3 independent transforms.
    The layout of the data is
    - all my own elements first (whether used or not)
    - followed by used-only remote elements sorted by remote processor.
    - followed by - for each transformation index - the set of local or
    remote elements with that transformation.
    The inputs for the constructor are
    - the set of untransformed local or remote indices in globalIndex
    numbering. These get rewritten to be indices into the layout of the data.
    - the set of transformed local or remote indices in globalIndexAndTransform
    encoding. These are labelPairs.

    Any distribute with transforms is now done as:
    1. exchange data with other processors and receive these into the
    slots for that processor
    2. for all transformations transform a subset of the data according
    to transformElements_[transformI] and store this starting from
    transformStart_[transformI]

    In the same way a reverse distribute will
    1. apply the inverse transform to the data starting at
    transformStart_[transformI] and copy the result back into the
    transformElements_[transformI]. These might be local or remote slots.
    2. the data in the remote slots will now be sent back to the correct
    location in the originating processor.

    E.g. a map to handle
    - mesh points on a mesh with
    - 1 cyclic so 3 permutations (+,-,none) will have layout
    - on e.g. processor 1 out of 2:

        +------+ <- transformStart[2]
        |      |
        |      | <- transform2 applied to data in local or remote slots
        |      |
        +------+ <- transformStart[1]
        |      |
        |      | <- transform1 applied to data in local or remote slots
        |      |
        +------+ <- transformStart[1]
        |      |
        |      | <- transform0 applied to data in local or remote slots
        |      |
        +------+ <- transformStart[0]
        |      |
        |      | <- data from proc2
        |      |
        +------+
        |      |
        |      | <- data from proc0
        |      |
        +------+ <- mesh.nPoints()
        |      |
        |      |
        |      |
        +------+ 0


SourceFiles
    mapDistribute.C
    mapDistributeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistribute_H
#define mapDistribute_H

#include "transformList.H"
#include "labelList.H"
#include "labelPair.H"
#include "Pstream.H"
#include "boolList.H"
#include "Map.H"
#include "vectorTensorTransform.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class globalIndex;
class PstreamBuffers;
class globalIndexAndTransform;

/*---------------------------------------------------------------------------*\
                           Class mapDistribute Declaration
\*---------------------------------------------------------------------------*/

class mapDistribute
{
    // Private data

        //- Size of reconstructed data
        label constructSize_;

        //- Maps from subsetted data back to original data
        labelListList subMap_;

        //- Maps from subsetted data to new reconstructed data
        labelListList constructMap_;

        // Optional transformation

            //- For every globalIndexAndTransform::transformPermutations
            //  gives the elements that need to be transformed
            labelListList transformElements_;

            //- Destination in constructMap for transformed elements
            labelList transformStart_;


        //- Schedule
        mutable autoPtr<List<labelPair> > schedulePtr_;


   // Private Member Functions

        static void checkReceivedSize
        (
            const label procI,
            const label expectedSize,
            const label receivedSize
        );

        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelList& elements,
            List<Map<label> >& compactMap
        ) const;
        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelListList& elements,
            List<Map<label> >& compactMap
        ) const;

        void exchangeAddressing
        (
            const int tag,
            const globalIndex& globalNumbering,
            labelList& elements,
            List<Map<label> >& compactMap,
            labelList& compactStart
        );
        void exchangeAddressing
        (
            const int tag,
            const globalIndex& globalNumbering,
            labelListList& elements,
            List<Map<label> >& compactMap,
            labelList& compactStart
        );


        //- Helper function: copy transformElements without transformation
        template<class T>
        void applyDummyTransforms(List<T>& field) const;

        template<class T, class TransformOp>
        void applyTransforms
        (
            const globalIndexAndTransform& globalTransforms,
            List<T>& field,
            const TransformOp& top
        ) const;

        //- Helper function: copy transformElements without transformation
        template<class T>
        void applyDummyInverseTransforms(List<T>& field) const;

        template<class T, class TransformOp>
        void applyInverseTransforms
        (
            const globalIndexAndTransform& globalTransforms,
            List<T>& field,
            const TransformOp& top
        ) const;


public:

    // Public classes

        //- Default transformation behaviour
        class transform
        {
        public:

            template<class Type>
            void operator()
            (
                const vectorTensorTransform& vt,
                const bool forward,
                List<Type>& fld
            ) const
            {
                const tensor T(forward ? vt.R() : vt.R().T());
                transformList(T, fld);
            }

            template<class Type>
            void operator()
            (
                const vectorTensorTransform& vt,
                const bool forward,
                List<List<Type> >& flds
            ) const
            {
                forAll(flds, i)
                {
                    operator()(vt, forward, flds[i]);
                }
            }

            //- Transform patch-based field
            template<class Type>
            void operator()(const coupledPolyPatch& cpp, UList<Type>& fld) const
            {
                if (!cpp.parallel())
                {
                    transformList(cpp.forwardT(), fld);
                }
            }

            //- Transform sparse field
            template<class Type, template<class> class Container>
            void operator()(const coupledPolyPatch& cpp, Container<Type>& map)
            const
            {
                if (!cpp.parallel())
                {
                    transformList(cpp.forwardT(), map);
                }
            }
        };

        //- Default transformation behaviour for position
        class transformPosition
        {
        public:

            void operator()
            (
                const vectorTensorTransform& vt,
                const bool forward,
                List<point>& fld
            ) const
            {
                pointField pfld(fld.xfer());
                if (forward)
                {
                    fld = vt.transformPosition(pfld);
                }
                else
                {
                    fld = vt.invTransformPosition(pfld);
                }
            }
            void operator()
            (
                const vectorTensorTransform& vt,
                const bool forward,
                List<List<point> >& flds
            ) const
            {
                forAll(flds, i)
                {
                    operator()(vt, forward, flds[i]);
                }
            }
            //- Transform patch-based field
            void operator()(const coupledPolyPatch& cpp, pointField& fld) const
            {
                cpp.transformPosition(fld);
            }
            template<template<class> class Container>
            void operator()(const coupledPolyPatch& cpp, Container<point>& map)
            const
            {
                Field<point> fld(map.size());
                label i = 0;
                forAllConstIter(typename Container<point>, map, iter)
                {
                    fld[i++] = iter();
                }
                cpp.transformPosition(fld);
                i = 0;
                forAllIter(typename Container<point>, map, iter)
                {
                    iter() = fld[i++];
                }
            }
        };



    // Declare name of the class and its debug switch
    ClassName("mapDistribute");


    // Constructors

        //- Construct null
        mapDistribute();

        //- Construct from components
        mapDistribute
        (
            const label constructSize,
            const Xfer<labelListList>& subMap,
            const Xfer<labelListList>& constructMap
        );

        //- Construct from components
        mapDistribute
        (
            const label constructSize,
            const Xfer<labelListList>& subMap,
            const Xfer<labelListList>& constructMap,
            const Xfer<labelListList>& transformElements,
            const Xfer<labelList>& transformStart
        );

        //- Construct from reverse addressing: per data item the send
        //  processor and the receive processor. (note: data is not stored
        //  sorted per processor so cannot use printLayout).
        mapDistribute
        (
            const labelList& sendProcs,
            const labelList& recvProcs
        );

        //- Construct from list of (possibly) remote elements in globalIndex
        //  numbering (or -1). Determines compact numbering (see above) and
        //  distribute map to get data into this ordering and renumbers the
        //  elements to be in compact numbering.
        mapDistribute
        (
            const globalIndex&,
            labelList& elements,
            List<Map<label> >& compactMap,
            const int tag = Pstream::msgType()
        );

        //- Special variant that works with the info sorted into bins
        //  according to local indices. E.g. think cellCells where
        //  cellCells[localCellI] is a list of global cells
        mapDistribute
        (
            const globalIndex&,
            labelListList& cellCells,
            List<Map<label> >& compactMap,
            const int tag = Pstream::msgType()
        );

        //- Construct from list of (possibly remote) untransformed elements
        //  in globalIndex numbering (or -1) and (possibly remote)
        //  transformded elements in globalIndexAndTransform numbering.
        //  Determines compact numbering (see above) and
        //  distribute map to get data into this ordering and renumbers the
        //  elements to be in compact numbering.
        mapDistribute
        (
            const globalIndex&,
            labelList& untransformedElements,
            const globalIndexAndTransform&,
            const labelPairList& transformedElements,
            labelList& transformedIndices,
            List<Map<label> >& compactMap,
            const int tag = Pstream::msgType()
        );

        //- As above but with ListLists.
        mapDistribute
        (
            const globalIndex&,
            labelListList& cellCells,
            const globalIndexAndTransform&,
            const List<labelPairList>& transformedElements,
            labelListList& transformedIndices,
            List<Map<label> >& compactMap,
            const int tag = Pstream::msgType()
        );

        //- Construct by transferring parameter content
        mapDistribute(const Xfer<mapDistribute>&);

        //- Construct copy
        mapDistribute(const mapDistribute&);


    // Member Functions

        // Access

            //- Constructed data size
            label constructSize() const
            {
                return constructSize_;
            }

            //- Constructed data size
            label& constructSize()
            {
                return constructSize_;
            }

            //- From subsetted data back to original data
            const labelListList& subMap() const
            {
                return subMap_;
            }

            //- From subsetted data back to original data
            labelListList& subMap()
            {
                return subMap_;
            }

            //- From subsetted data to new reconstructed data
            const labelListList& constructMap() const
            {
                return constructMap_;
            }

            //- From subsetted data to new reconstructed data
            labelListList& constructMap()
            {
                return constructMap_;
            }

            //- For every globalIndexAndTransform::transformPermutations
            //  gives the elements that need to be transformed
            const labelListList& transformElements() const
            {
                return transformElements_;
            }

            //- Destination in constructMap for transformed elements
            const labelList& transformStart() const
            {
                return transformStart_;
            }

            //- Find transform from transformElements
            label whichTransform(const label index) const;

            //- Calculate a schedule. See above.
            static List<labelPair> schedule
            (
                const labelListList& subMap,
                const labelListList& constructMap,
                const int tag
            );

            //- Return a schedule. Demand driven. See above.
            const List<labelPair>& schedule() const;


        // Other

            //- Transfer the contents of the argument and annul the argument.
            void transfer(mapDistribute&);

            //- Transfer contents to the Xfer container
            Xfer<mapDistribute> xfer();

            //- Helper for construct from globalIndex. Renumbers element
            //  (in globalIndex numbering) into compact indices.
            static label renumber
            (
                const globalIndex&,
                const List<Map<label> >& compactMap,
                const label globalElement
            );

            //- Compact maps. Gets per field a bool whether it is used (locally)
            //  and works out itself what this side and sender side can remove
            //  from maps.
            void compact
            (
                const boolList& elemIsUsed,
                const int tag = UPstream::msgType()
            );

            //- Distribute data. Note:schedule only used for Pstream::scheduled
            //  for now, all others just use send-to-all, receive-from-all.
            template<class T>
            static void distribute
            (
                const Pstream::commsTypes commsType,
                const List<labelPair>& schedule,
                const label constructSize,
                const labelListList& subMap,
                const labelListList& constructMap,
                List<T>&,
                const int tag = UPstream::msgType()
            );

            //- Distribute data. If multiple processors writing to same
            //  position adds contributions using cop.
            template<class T, class CombineOp>
            static void distribute
            (
                const Pstream::commsTypes commsType,
                const List<labelPair>& schedule,
                const label constructSize,
                const labelListList& subMap,
                const labelListList& constructMap,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const int tag = UPstream::msgType()
            );

            //- Distribute data using default commsType.
            template<class T>
            void distribute
            (
                DynamicList<T>& fld,
                const bool dummyTransform = true,
                const int tag = UPstream::msgType()
            ) const;

            //- Distribute data using default commsType.
            template<class T>
            void distribute
            (
                List<T>& fld,
                const bool dummyTransform = true,
                const int tag = UPstream::msgType()
            ) const;

            //- Same but with transforms
            template<class T, class TransformOp>
            void distribute
            (
                const globalIndexAndTransform&,
                List<T>& fld,
                const TransformOp& top,
                const int tag = UPstream::msgType()
            ) const;

            //- Reverse distribute data using default commsType.
            template<class T>
            void reverseDistribute
            (
                const label constructSize,
                List<T>&,
                const bool dummyTransform = true,
                const int tag = UPstream::msgType()
            ) const;

            //- Same but with transforms
            template<class T, class TransformOp>
            void reverseDistribute
            (
                const globalIndexAndTransform&,
                const label constructSize,
                List<T>& fld,
                const TransformOp& top,
                const int tag = UPstream::msgType()
            ) const;

            //- Reverse distribute data using default commsType.
            //  Since constructSize might be larger than supplied size supply
            //  a nullValue
            template<class T>
            void reverseDistribute
            (
                const label constructSize,
                const T& nullValue,
                List<T>& fld,
                const bool dummyTransform = true,
                const int tag = UPstream::msgType()
            ) const;

            //- Same but with transforms
            template<class T, class TransformOp>
            void reverseDistribute
            (
                const globalIndexAndTransform&,
                const label constructSize,
                const T& nullValue,
                List<T>& fld,
                const TransformOp& top,
                const int tag = UPstream::msgType()
            ) const;

            //- Do all sends using PstreamBuffers
            template<class T>
            void send(PstreamBuffers&, const List<T>&) const;
            //- Do all receives using PstreamBuffers
            template<class T>
            void receive(PstreamBuffers&, List<T>&) const;

            //- Debug: print layout. Can only be used on maps with sorted
            //  storage (local data first, then non-local data)
            void printLayout(Ostream& os) const;

            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistribute::updateMesh(const mapPolyMesh&)"
                );
            }

    // Member Operators

        void operator=(const mapDistribute&);


    // IOstream operators

        //- Read dictionary from Istream
        friend Istream& operator>>(Istream&, mapDistribute&);

        //- Write dictionary to Ostream
        friend Ostream& operator<<(Ostream&, const mapDistribute&);

};


template<>
void mapDistribute::transform::operator()
(
    const vectorTensorTransform&,
    const bool,
    List<label>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    UList<label>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    Map<label>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<label>&
) const;

template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    UList<scalar>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const vectorTensorTransform&,
    const bool,
    List<scalar>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    Map<scalar>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<scalar>&
) const;

template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch& cpp,
    UList<bool>& fld
) const;
template<>
void mapDistribute::transform::operator()
(
    const vectorTensorTransform&,
    const bool,
    List<bool>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    Map<bool>&
) const;
template<>
void mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<bool>&
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mapDistributeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
