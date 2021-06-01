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
    Foam::PackedList

Description
    A dynamically allocatable list of packed unsigned integers.

    The list resizing is similar to DynamicList, thus the methods clear()
    and setSize() behave like their DynamicList counterparts and the methods
    reserve() and setCapacity() can be used to influence the allocation.

    The number of bits per item is specified by the template parameter nBits.

Note
    In a const context, the '[]' operator simply returns the stored value,
    with out-of-range elements returned as zero.
    In a non-const context, the '[]' operator returns an iteratorBase, which
    might not have a valid reference for out-of-range elements.
    The iteratorBase class handles the assignment of new values.

    Using the iteratorBase as a proxy allows assignment of values
    between list elements. Thus the following bit of code works as expected:
    \code
        list[1] = list[5];      // value assignment, not iterator position
        list[2] = list[5] = 4;  // propagates value
        list[1] = list[5] = list[6];  // propagates value
    \endcode

    Using get() or the '[]' operator are similarly fast. Looping and reading
    via an iterator is approx. 15% slower, but can be more flexible.

    Using the set() operator (and the '[]' operator) are marginally slower
    (approx. 5%) than using an iterator, but the set() method has the
    advantage of also returning a bool if the value changed.  This can be
    useful for branching on changed values.

    \code
        list[5] = 4;
        changed = list.set(5, 8);
        if (changed) ...
    \endcode

    The lazy evaluation used means that reading an out-of-range element
    returns zero, but does not affect the list size.  Even in a non-const
    context, only the assignment itself causes the element to be created.
    For example,
    \code
        list.resize(4);
        Info<< list[10] << "\n";  // print zero, but doesn't adjust list
        list[8] = 1;
    \endcode

    Also note that all unused internal storage elements are guaranteed to
    always be bit-wise zero. This property must not be violated by any
    inheriting classes.

    In addition to the normal output format, PackedList also supports a
    compact ASCII format that may be convenient for user input in some
    situations. The general format is a group of index/value pairs:
    \verbatim
        { (index1 value1) (index2 value2) (index3 value3) }
    \endverbatim
    The bool specialisation just uses the indices corresponding to
    non-zero entries instead of a index/value pair:
    \verbatim
        { index1 index2 index3 }
    \endverbatim
    In both cases, the supplied indices can be randomly ordered.

See also
    Foam::DynamicList

SourceFiles
    PackedListI.H
    PackedList.C

\*---------------------------------------------------------------------------*/

#ifndef PackedList_H
#define PackedList_H

#include "labelList.H"
#include "UIndirectList.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<unsigned nBits> class PackedList;

template<unsigned nBits>
void writeEntry(Ostream& os, const PackedList<nBits>&);

template<unsigned nBits>
Istream& operator>>(Istream&, PackedList<nBits>&);

template<unsigned nBits>
Ostream& operator<<(Ostream&, const PackedList<nBits>&);


/*---------------------------------------------------------------------------*\
                       Class PackedListCore Declaration
\*---------------------------------------------------------------------------*/

//- Template-invariant bits for PackedList
struct PackedListCore
{
    //- Construct null
    PackedListCore()
    {}

    //- Define template name and debug
    ClassName("PackedList");
};


/*---------------------------------------------------------------------------*\
                         Class PackedList Declaration
\*---------------------------------------------------------------------------*/

template<unsigned nBits=1>
class PackedList
:
    public PackedListCore,
    private List<unsigned int>
{
protected:

    typedef unsigned int      StorageType;
    typedef List<StorageType> StorageList;

    // Protected Member Functions

        //- Calculate the list length when packed
        inline static label packedLength(const label);

        //- Read a list entry (allows for specialisation)
        inline static unsigned int readValue(Istream&);

        //- Read an index/value pair and set accordingly.
        //  For bool specialisation, read a single index value
        inline void setPair(Istream&);


private:

    //- nBits must be positive (non-zero) and fit within the storage.
    //  For efficiency, however, require packing at least 2 items otherwise
    //  it is more efficient to use a normal list.
    //  Thus max nBits is 1/2 of the base storage size.
    //  For simplicity, assume 8-bit bytes in the assert.
    static_assert
    (
        nBits && nBits <= (sizeof(StorageType) << 2),
        "nBits must be positive (non-zero) and fit within the storage"
    );

    // Private Data

        //- Number of nBits entries
        label size_;


public:

    // Public data

        //- The max. number of bits that can be templated.
        //  Might someday be useful for a template assert.
        inline static unsigned int max_bits();

        //- The max. value for an entry, which simultaneously the bit-mask
        //  eg, ((1 << 2) - 1) yields 0b0011
        inline static unsigned int max_value();

        //- The number of entries per packed storage element
        inline static unsigned int packing();

        //- Masking for all bits below the offset
        inline static unsigned int maskLower(unsigned offset);


    // Forward declaration of iterators

        class iteratorBase;
        class iterator;
        class const_iterator;


    // Constructors

        //- Null constructor
        inline PackedList();

        //- Construct with given size, initialises list to 0
        explicit inline PackedList(const label size);

        //- Construct with given size and value for all elements
        inline PackedList(const label size, const unsigned val);

        //- Construct from Istream
        inline PackedList(Istream&);

        //- Copy constructor
        inline PackedList(const PackedList<nBits>&);

        //- Move constructor
        inline PackedList(PackedList<nBits>&&);

        //- Construct from a list of labels
        explicit inline PackedList(const labelUList&);

        //- Construct from an indirect list of labels
        explicit inline PackedList(const UIndirectList<label>&);

        //- Clone
        inline autoPtr<PackedList<nBits>> clone() const;


    // Member Functions

        // Access

            //- The number of elements that can be stored before reallocating
            inline label capacity() const;

            //- Number of entries.
            inline label size() const;

            //- Return true if the list is empty (ie, size() is zero).
            inline bool empty() const;

            //- Get value at index I.
            //  Never auto-vivify entries.
            inline unsigned int get(const label) const;

            //- Set value at index I. Return true if value changed.
            //  Does auto-vivify for non-existent entries.
            //  Default value set is the max_value.
            inline bool set(const label, const unsigned int val = ~0u);

            //- Unset the entry at index I. Return true if value changed.
            //  Never auto-vivify entries.
            inline bool unset(const label);

            //- Return the underlying packed storage
            //  Manipulate with utmost caution
            inline List<unsigned int>& storage();

            //- Return the underlying packed storage
            inline const List<unsigned int>& storage() const;

            //- The list length when packed
            inline label packedLength() const;

            //- Return the binary size in number of characters
            //  used in the underlying storage
            inline std::streamsize byteSize() const;

            //- Count number of bits set, O(log(n))
            //  Uses the Hamming weight (population count) method
            //  http://en.wikipedia.org/wiki/Hamming_weight
            unsigned int count() const;

            //- Return the values as a list of labels
            labelList values() const;

            //- Print bit patterns, optionally output unused elements
            //
            // addressable bits:
            //   on: '1', off: '-'
            //
            // non-addressable bits:
            //   on: '!', off: '.'
            //
            Ostream& printBits(Ostream&, const bool fullOutput=false) const;

            //- Print information and bit patterns (with printBits)
            Ostream& printInfo(Ostream&, const bool fullOutput=false) const;


        // Edit

            //- Trim any trailing zero elements
            bool trim();

            //- Invert the bits in the addressable region
            void flip();

            //- Clear all bits
            inline void reset();

            //- Alter the size of the underlying storage.
            //  The addressed size will be truncated if needed to fit, but will
            //  remain otherwise untouched.
            inline void setCapacity(const label);

            //- Reset addressable list size, does not shrink the allocated size.
            //  Optionally specify a value for new elements.
            inline void resize(const label, const unsigned int& val = 0u);

            //- Alias for resize()
            inline void setSize(const label, const unsigned int& val = 0u);

            //- Reserve allocation space for at least this size.
            //  Never shrinks the allocated size.
            //  The list size is adjusted as per DynamicList with
            //  SizeInc=0, SizeMult=2, SizeDiv=1
            inline void reserve(const label);

            //- Clear the list, i.e. set addressable size to zero.
            //  Does not adjust the underlying storage
            inline void clear();

            //- Clear the list and delete storage.
            inline void clearStorage();

            //- Shrink the allocated space to what is actually used.
            inline void shrink();

            //- Transfer the contents of the argument list into this list
            //  and annul the argument list.
            inline void transfer(PackedList<nBits>&);


        // IO

            //- Clear list and read from stream
            Istream& read(Istream&);

            //- Write, optionally with indexedOutput
            //
            // The indexed output may be convenient in some situations.
            // The general format is a group of index/value pairs:
            //  \verbatim
            //      { (index1 value1) (index2 value2) (index3 value3) }
            // \endverbatim
            // The bool specialisation just uses the indices corresponding to
            // non-zero entries instead of a index/value pair:
            // \verbatim
            //     { index1 index2 index3 }
            // \endverbatim
            //
            // Note the indexed output is only supported for ASCII streams.
            Ostream& write
            (
                Ostream&,
                const bool indexedOutput=false
            ) const;


    // Member Operators

            //- Append a value at the end of the list
            inline PackedList<nBits>& append(const unsigned int val);

            //- Remove and return the last element
            inline unsigned int remove();

            //- Get value at index I
            //  Never auto-vivify entries.
            inline unsigned int operator[](const label) const;

            //- Set value at index I.
            //  Returns iterator to perform the actual operation.
            //  Does not auto-vivify entries, but will when assigned to.
            inline iteratorBase operator[](const label);

            //- Assignment of all entries to the given value. Takes linear time.
            inline void operator=(const unsigned int val);

            //- Assignment operator
            void operator=(const PackedList<nBits>&);

            //- Move assignment operator
            void operator=(PackedList<nBits>&&);

            //- Assignment operator
            void operator=(const labelUList&);

            //- Assignment operator
            void operator=(const UIndirectList<label>&);


    // Iterators and helpers

        //- The iterator base for PackedList
        //  Note: data and functions are protected, to allow reuse by iterator
        //  and prevent most external usage.
        class iteratorBase
        {
            friend class PackedList;

        protected:

            // Protected Data

                //- Pointer to original list
                //  This also lets us use the default bitwise copy/assignment
                PackedList* list_;

                //- Element index
                label index_;


            // Protected Member Functions

                //- Get value as unsigned, no range-checking
                inline unsigned int get() const;

                //- Set value, returning true if changed, no range-checking
                inline bool set(unsigned int);


            // Constructors

                //- Construct null
                inline iteratorBase();

                //- Construct from base list and position index
                inline iteratorBase(const PackedList*, const label);

                //- Default copy constructor
                iteratorBase(const iteratorBase&) = default;


        public:

            // Member Functions

                //- Return the element index corresponding to the iterator
                inline label key() const;

                //- Write index/value for a non-zero entry
                //  The bool specialisation writes the index only
                inline bool writeIfSet(Ostream&) const;

            // Member Operators

                //- Compare values (not positions)
                inline bool operator==(const iteratorBase&) const;
                inline bool operator!=(const iteratorBase&) const;

                //- Assign value, not position.
                //  This allows packed[0] = packed[3] for assigning values
                inline void operator=(const iteratorBase&);

                //- Assign value.
                //  A non-existent entry will be auto-vivified.
                inline void operator=(const unsigned int val);

                //- Conversion operator
                //  Never auto-vivify entries.
                inline operator unsigned int () const;

            //- Print information and values
            Ostream& printInfo(Ostream&) const;
        };


        //- The iterator class used for PackedList
        class iterator
        :
            public iteratorBase
        {

            //- Disallow default bitwise copy construction
            //  This would violate const-ness!
            iterator(const const_iterator&);

            //- Disallow assignment from const_iterator
            //  This would violate const-ness!
            void operator=(const const_iterator&);


        public:

            // Constructors

                //- Construct null
                inline iterator();

                //- Construct from iterator base, eg iter(packedlist[i])
                //  but also  "iterator iter = packedlist[i];"
                //  An out-of-range iterator is assigned end()
                inline iterator(const iteratorBase&);

                //- Construct from base list and position index
                inline iterator(const PackedList*, const label);


            // Member Operators

                //- Compare positions (not values)
                inline bool operator==(const iteratorBase&) const;
                inline bool operator!=(const iteratorBase&) const;

                //- Assign from iteratorBase, eg iter = packedlist[i]
                //  An out-of-range iterator is assigned end()
                inline void operator=(const iteratorBase&);

                //- Return value
                inline unsigned int operator*() const;

                //- Return value
                inline unsigned int operator()() const;

                //- Return iteratorBase for assigning values
                inline iteratorBase& operator*();

                //- Return iteratorBase for assigning values
                inline iteratorBase& operator()();

                inline iterator& operator++();
                inline iterator operator++(int);

                inline iterator& operator--();
                inline iterator operator--(int);

        };

        //- Iterator set to the beginning of the PackedList
        inline iterator begin();

        //- Iterator set to beyond the end of the PackedList
        inline iterator end();


        //- The const_iterator for PackedList
        class const_iterator
        :
            public iteratorBase
        {
        public:

            // Constructors

                //- Construct null
                inline const_iterator();

                //- Construct from iterator base, eg iter(packedlist[i])
                //  but also  "const_iterator iter = packedlist[i];"
                //  An out-of-range iterator is assigned cend()
                inline const_iterator(const iteratorBase&);

                //- Construct from base list and position index
                inline const_iterator(const PackedList*, const label);

                //- Construct from iterator
                inline const_iterator(const iterator&);


            // Member Operators

                //- Compare positions (not values)
                inline bool operator==(const iteratorBase&) const;
                inline bool operator!=(const iteratorBase&) const;

                //- Assign from iteratorBase or derived
                //  eg, iter = packedlist[i] or even iter = list.begin()
                inline void operator=(const iteratorBase&);

                //- Return referenced value directly
                inline unsigned int operator*() const;

                //- Return referenced value directly
                inline unsigned int operator()() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);

                inline const_iterator& operator--();
                inline const_iterator operator--(int);
        };


        //- const_iterator set to the beginning of the PackedList
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the PackedList
        inline const_iterator cend() const;

        //- const_iterator set to the beginning of the PackedList
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the PackedList
        inline const_iterator end() const;


    // IOstream Operators

        friend Istream& operator>> <nBits>
        (
            Istream&,
            PackedList<nBits>&
        );

        friend Ostream& operator<< <nBits>
        (
            Ostream&,
            const PackedList<nBits>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PackedListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PackedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
