/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::PackedBoolList

Description
    A bit-packed bool list.

    In addition to the obvious memory advantage over using a
    List\<bool\>, this class also provides a number of bit-like
    operations.

SourceFiles
    PackedBoolListI.H
    PackedBoolList.C

See also
    Foam::PackedList

\*---------------------------------------------------------------------------*/

#ifndef PackedBoolList_H
#define PackedBoolList_H

#include "PackedList.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class PackedBoolList;

//- \typedef A List of PackedBoolList
typedef List<PackedBoolList> PackedBoolListList;

/*---------------------------------------------------------------------------*\
                       Class PackedBoolList Declaration
\*---------------------------------------------------------------------------*/

class PackedBoolList
:
    public PackedList<1>
{
    // Private Member Functions

        //- Preparation, resizing before a bitor operation
        //  returns true if the later result needs trimming
        bool bitorPrepare(const PackedList<1>& lst, label& maxPackLen);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        template<class LabelListType>
        label setIndices(const LabelListType& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        template<class LabelListType>
        label unsetIndices(const LabelListType& indices);

        //- Subset with the listed indices. Return number of elements subsetted.
        template<class LabelListType>
        label subsetIndices(const LabelListType& indices);


public:

    // Constructors

        //- Construct null
        inline PackedBoolList();

        //- Construct from Istream
        PackedBoolList(Istream&);

        //- Construct with given size, initializes list to 0
        explicit inline PackedBoolList(const label size);

        //- Construct with given size and value for all elements
        inline PackedBoolList(const label size, const bool val);

        //- Copy constructor
        inline PackedBoolList(const PackedBoolList&);

        //- Copy constructor
        explicit inline PackedBoolList(const PackedList<1>&);

        //- Construct by transferring the parameter contents
        inline PackedBoolList(const Xfer<PackedBoolList>&);

        //- Construct by transferring the parameter contents
        inline PackedBoolList(const Xfer<PackedList<1>>&);

        //- Construct from a list of bools
        explicit inline PackedBoolList(const Foam::UList<bool>&);

        //- Construct from a list of labels
        //  using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const labelUList& indices);

        //- Construct from a list of labels
        //  using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const UIndirectList<label>& indices);

        //- Clone
        inline autoPtr<PackedBoolList> clone() const;


    // Member Functions

        // Access

            using PackedList<1>::set;
            using PackedList<1>::unset;

            //- Set specified bits.
            void set(const PackedList<1>&);

            //- Set the listed indices. Return number of elements changed.
            //  Does auto-vivify for non-existent entries.
            label set(const labelUList& indices);

            //- Set the listed indices. Return number of elements changed.
            //  Does auto-vivify for non-existent entries.
            label set(const UIndirectList<label>& indices);

            //- Unset specified bits.
            void unset(const PackedList<1>&);

            //- Unset the listed indices. Return number of elements changed.
            //  Never auto-vivify entries.
            label unset(const labelUList& indices);

            //- Unset the listed indices. Return number of elements changed.
            //  Never auto-vivify entries.
            label unset(const UIndirectList<label>& indices);

            //- Subset with the specified list.
            void subset(const PackedList<1>&);

            //- Subset with the listed indices.
            //  Return number of elements subsetted.
            label subset(const labelUList& indices);

            //- Subset with the listed indices.
            //  Return number of elements subsetted.
            label subset(const UIndirectList<label>& indices);


            //- Return indices of the used (true) elements as a list of labels
            Xfer<labelList> used() const;


        // Edit

            //- Transfer the contents of the argument list into this list
            //  and annul the argument list.
            inline void transfer(PackedBoolList&);

            //- Transfer the contents of the argument list into this list
            //  and annul the argument list.
            inline void transfer(PackedList<1>&);

            //- Transfer contents to the Xfer container
            inline Xfer<PackedBoolList> xfer();


    // Member Operators

            //- Assignment of all entries to the given value.
            inline void operator=(const bool val);

            //- Assignment operator.
            inline void operator=(const PackedBoolList&);

            //- Assignment operator.
            inline void operator=(const PackedList<1>&);

            //- Assignment operator.
            void operator=(const Foam::UList<bool>&);

            //- Assignment operator,
            //  using the labels as indices to indicate which bits are set
            inline void operator=(const labelUList& indices);

            //- Assignment operator,
            //  using the labels as indices to indicate which bits are set
            inline void operator=(const UIndirectList<label>&);

            //- Complement operator
            inline PackedBoolList operator~() const;

            //- And operator (lists may be dissimilar sizes)
            inline PackedBoolList& operator&=(const PackedList<1>&);

            //- And operator (lists may be dissimilar sizes)
            //  using the labels as indices to indicate which bits are set
            inline PackedBoolList& operator&=(const labelUList& indices);

            //- And operator (lists may be dissimilar sizes)
            //  using the labels as indices to indicate which bits are set
            inline PackedBoolList& operator&=(const UIndirectList<label>&);

            //- Xor operator (lists may be dissimilar sizes)
            //  Retains unique entries
            PackedBoolList& operator^=(const PackedList<1>&);

            //- Or operator (lists may be dissimilar sizes)
            inline PackedBoolList& operator|=(const PackedList<1>&);

            //- Or operator (lists may be dissimilar sizes),
            //  using the labels as indices to indicate which bits are set
            inline PackedBoolList& operator|=(const labelUList& indices);

            //- Or operator (lists may be dissimilar sizes),
            //  using the labels as indices to indicate which bits are set
            inline PackedBoolList& operator|=(const UIndirectList<label>&);


            //- Add entries to this list, synonymous with the or operator
            inline PackedBoolList& operator+=(const PackedList<1>&);

            //- Add entries to this list, synonymous with the or operator
            inline PackedBoolList& operator+=(const labelUList& indices);

            //- Add entries to this list, synonymous with the or operator
            inline PackedBoolList& operator+=(const UIndirectList<label>&);

            //- Remove entries from this list - unset the specified bits
            inline PackedBoolList& operator-=(const PackedList<1>&);

            //- Remove entries from this list - unset the specified bits
            inline PackedBoolList& operator-=(const labelUList& indices);

            //- Remove entries from this list - unset the specified bits
            inline PackedBoolList& operator-=(const UIndirectList<label>&);
};


// Global Operators

//- Intersect lists - the result is trimmed to the smallest intersecting size
PackedBoolList operator&
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


//- Combine to form a unique list (xor)
//  The result is trimmed to the smallest intersecting size
PackedBoolList operator^
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


//- Combine lists
PackedBoolList operator|
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PackedBoolListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
