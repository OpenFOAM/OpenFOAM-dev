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
    Foam::List

Description
    A 1D array of objects of type \<T\>, where the size of the vector
    is known and used for subscript bounds checking, etc.

    Storage is allocated on free-store during construction.

SourceFiles
    List.C
    ListI.H
    ListIO.C

\*---------------------------------------------------------------------------*/

#ifndef List_H
#define List_H

#include "UList.H"
#include "autoPtr.H"
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

using std::move;

template<class T>
inline T clone(const T& t)
{
    return move(T(t));
}

// Forward declaration of classes

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class T> class List;

template<class T> Istream& operator>>(Istream&, List<T>&);

template<class T, unsigned Size> class FixedList;
template<class T> class PtrList;

class SLListBase;
template<class LListBase, class T> class LList;
template<class T>
using SLList = LList<SLListBase, T>;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicList;

template<class T> class SortableList;
template<class T> class IndirectList;
template<class T> class UIndirectList;
template<class T> class BiIndirectList;

typedef UList<label> unallocLabelList;

/*---------------------------------------------------------------------------*\
                            Class List Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class List
:
    public UList<T>
{
    // Private Member Functions

        //- Allocate list storage
        inline void alloc();

        //- Reallocate list storage to the given size
        inline void reAlloc(const label s);

        //- Copy list of given type
        template<class List2>
        inline void copyList(const List2&);

        //- Allocate storage and copy list of given type
        template<class List2>
        inline void allocCopyList(const List2&);

        //- Construct given start and end iterators and number of elements
        template<class InputIterator>
        inline List(InputIterator first, InputIterator last, const label s);


protected:

    //- Override size to be inconsistent with allocated storage.
    //  Use with care
    inline void size(const label);


public:

    // Static Member Functions

        //- Return a null List
        inline static const List<T>& null();


    // Constructors

        //- Null constructor
        inline List();

        //- Construct with given size
        explicit List(const label);

        //- Construct with given size and value for all elements
        List(const label, const T&);

        //- Construct with given size initialising all elements to zero
        List(const label, const zero);

        //- Copy constructor
        List(const List<T>&);

        //- Copy constructor from list containing another type
        template<class T2>
        explicit List(const List<T2>&);

        //- Move constructor
        List(List<T>&&);

        //- Construct as copy or re-use as specified
        List(List<T>&, bool reuse);

        //- Construct as subset
        List(const UList<T>&, const labelUList& mapAddressing);

        //- Construct given start and end iterators
        template<class InputIterator>
        List(InputIterator first, InputIterator last);

        //- Construct as copy of FixedList<T, Size>
        template<unsigned Size>
        explicit List(const FixedList<T, Size>&);

        //- Construct as copy of PtrList<T>
        explicit List(const PtrList<T>&);

        //- Construct as copy of SLList<T>
        explicit List(const SLList<T>&);

        //- Construct as copy of UIndirectList<T>
        explicit List(const UIndirectList<T>&);

        //- Construct as copy of BiIndirectList<T>
        explicit List(const BiIndirectList<T>&);

        //- Construct from an initialiser list
        List(std::initializer_list<T>);

        //- Construct from Istream
        List(Istream&);

        //- Clone
        inline autoPtr<List<T>> clone() const;


    //- Destructor
    ~List();


    // Related types

        //- Declare type of subList
        typedef SubList<T> subList;


    // Member Functions

        //- Return the number of elements in the UList
        inline label size() const;


        // Edit

            //- Alias for setSize(const label)
            inline void resize(const label);

            //- Alias for setSize(const label, const T&)
            inline void resize(const label, const T&);

            //- Reset size of List
            void setSize(const label);

            //- Reset size of List and value for new elements
            void setSize(const label, const T&);

            //- Clear the list, i.e. set size to zero
            inline void clear();

            //- Append an element at the end of the list
            inline void append(const T&);

            //- Append a List at the end of this list
            inline void append(const UList<T>&);

            //- Append a UIndirectList at the end of this list
            inline void append(const UIndirectList<T>&);

            //- Transfer the contents of the argument List into this list
            //  and annul the argument list
            void transfer(List<T>&);

            //- Transfer the contents of the argument List into this list
            //  and annul the argument list
            template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
            void transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>&);

            //- Transfer the contents of the argument List into this list
            //  and annul the argument list
            void transfer(SortableList<T>&);

            //- Return subscript-checked element of UList
            inline T& newElmt(const label);


        //- Disallow implicit shallowCopy
        void shallowCopy(const UList<T>&) = delete;


    // Member Operators

        //- Assignment to UList operator. Takes linear time
        void operator=(const UList<T>&);

        //- Assignment operator. Takes linear time
        void operator=(const List<T>&);

        //- Move assignment operator
        void operator=(List<T>&&);

        //- Assignment to SLList operator. Takes linear time
        void operator=(const SLList<T>&);

        //- Assignment to UIndirectList operator. Takes linear time
        void operator=(const UIndirectList<T>&);

        //- Assignment to BiIndirectList operator. Takes linear time
        void operator=(const BiIndirectList<T>&);

        //- Assignment to an initialiser list
        void operator=(std::initializer_list<T>);

        //- Assignment of all entries to the given value
        inline void operator=(const T&);

        //- Assignment of all entries to zero
        inline void operator=(const zero);


    // Istream operator

        //- Read List from Istream, discarding contents of existing List
        friend Istream& operator>> <T>
        (Istream&, List<T>&);
};


//- Read a bracket-delimited list, or handle a single value as list of size 1
//  For example,
//  \code
//      wList = readList<word>(IStringStream("(patch1 patch2 patch3)")());
//      wList = readList<word>(IStringStream("patch0")());
//  \endcode
//  Mostly useful for handling command-line arguments
template<class T>
List<T> readList(Istream&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "List.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
