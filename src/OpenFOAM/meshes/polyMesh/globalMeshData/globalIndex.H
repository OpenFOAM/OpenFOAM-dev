/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::globalIndex

Description
    Calculates a unique integer (label so might not have enough room - 2G max)
    for processor + local index. E.g.

    globalIndex globalFaces(mesh.nFaces());
    label globalFaceI = globalFaces.toGlobal(faceI);


SourceFiles
    globalIndexI.H
    globalIndex.C
    globalIndexTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef globalIndex_H
#define globalIndex_H

#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class globalIndex;

Istream& operator>>(Istream& is, globalIndex& gi);
Ostream& operator<<(Ostream& os, const globalIndex& gi);


/*---------------------------------------------------------------------------*\
                           Class globalIndex Declaration
\*---------------------------------------------------------------------------*/

class globalIndex
{
    // Private data

        //- Start of procI. Size is nProcs()+1. (so like CompactListList)
        labelList offsets_;


public:

    // Constructors

        //- Construct null
        inline globalIndex();

        //- Construct from local max size. Does communication with default
        //  communicator and message tag.
        globalIndex(const label localSize);

        //- Construct from local max size. Does communication with given
        //  communicator and message tag
        globalIndex
        (
            const label localSize,
            const int tag,
            const label comm,
            const bool parallel     // use parallel comms
        );

        //- Construct from components
        globalIndex(const labelList& offsets);

        //- Construct from components
        inline globalIndex(const Xfer<labelList>& offsets);

        //- Construct from Istream
        globalIndex(Istream& is);


    // Member Functions

        // Edit

            //- Change after construction
            inline labelList& offsets();


        // Queries relating to my processor (using world communicator)

            //- my local size
            inline label localSize() const;

            //- From local to global
            inline label toGlobal(const label i) const;

            //- Is on local processor
            inline bool isLocal(const label i) const;

            //- From global to local on current processor.
            //  FatalError if not on local processor.
            inline label toLocal(const label i) const;


        // Global queries

            //- Global sum of localSizes
            inline label size() const;

            //- Size of procI data
            inline label localSize(const label procI) const;

            //- From local to global on procI
            inline label toGlobal(const label procI, const label i) const;

            //- Is on processor procI
            inline bool isLocal(const label procI, const label i) const;

            //- From global to local on procI
            inline label toLocal(const label procI, const label i) const;

            //- Which processor does global come from? Binary search.
            inline label whichProcID(const label i) const;

            //- Start of procI data
            inline label offset(const label procI) const;


        // Other

            //- Collect data in processor order on master (== procIDs[0]).
            //  Needs offsets only on master.
            template<class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const labelList& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType=Pstream::nonBlocking
            );

            //- Collect data in processor order on master (== procIDs[0]).
            //  Needs offsets only on master.
            template<class Type>
            void gather
            (
                const label comm,
                const labelList& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType=Pstream::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, allFld, tag, commsType);
            }

            //- Inplace collect data in processor order on master
            //  (== procIDs[0]). Needs offsets only on master.
            template<class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const labelList& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType=Pstream::nonBlocking
            );

            //- Inplace collect data in processor order on master
            //  (== procIDs[0]). Needs offsets only on master.
            template<class Type>
            void gather
            (
                const label comm,
                const labelList& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType=Pstream::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, tag, commsType);
            }

            //- Distribute data in processor order. Requires fld to be sized!
            template<class Type>
            static void scatter
            (
                const labelUList& offsets,
                const label comm,
                const labelList& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType=Pstream::nonBlocking
            );

            //- Distribute data in processor order. Requires fld to be sized!
            template<class Type>
            void scatter
            (
                const label comm,
                const labelList& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType=Pstream::nonBlocking
            ) const
            {
                scatter(offsets_, comm, procIDs, allFld, fld, tag, commsType);
            }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, globalIndex& gi);
        friend Ostream& operator<<(Ostream& os, const globalIndex& gi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "globalIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalIndexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
