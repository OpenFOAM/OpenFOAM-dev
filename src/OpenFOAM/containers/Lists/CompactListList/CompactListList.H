/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::CompactListList

Description
    A packed storage unstructured matrix of objects of type \<T\>
    using an offset table for access.

    The offset table is the size of the number of rows+1
    whose elements are the
    accumulated sizes of the rows, i.e.
      - offset[i] gives the index of first element of row i
      - offset[i+1] - offset[i] is the number of elements in row i

    Storage is allocated on free-store during construction.

    As a special case a null-constructed CompactListList has an empty
    offsets_ (instead of size 1).

SourceFiles
    CompactListList.C
    CompactListListI.H
    CompactListListIO.C

\*---------------------------------------------------------------------------*/

#ifndef CompactListList_H
#define CompactListList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, class Container> class CompactListList;

template<class T, class Container> Istream& operator>>
(
    Istream&,
    CompactListList<T, Container>&
);
template<class T, class Container> Ostream& operator<<
(
    Ostream&,
    const CompactListList<T, Container>&
);


/*---------------------------------------------------------------------------*\
                       Class CompactListList Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Container = List<T>>
class CompactListList
{
    // Private Data

        label size_;

        //- Offset table
        List<label> offsets_;

        //- Packed matrix of data
        List<T> m_;


public:

    // Static Member Functions

        //- Return a null CompactListList
        inline static const CompactListList<T, Container>& null();


    // Constructors

        //- Null constructor.
        inline CompactListList();

        //- Move constructor
        CompactListList(CompactListList<T, Container>&&);

        //- Construct by converting given List<List<T>>
        explicit CompactListList(const List<Container>&);

        //- Construct given size of offset table (number of rows)
        //  and number of data.
        inline CompactListList(const label mRows, const label nData);

        //- Construct given size of offset table (number of rows),
        //  the number of data and a value for all elements.
        inline CompactListList(const label mRows, const label nData, const T&);

        //- Construct given list of row-sizes.
        explicit CompactListList(const labelUList& rowSizes);

        //- Construct given list of row-sizes
        CompactListList(const labelUList& rowSizes, const T&);

        //- Construct as copy or re-use as specified.
        CompactListList(CompactListList<T, Container>&, bool reuse);

        //- Construct from Istream.
        CompactListList(Istream&);

        //- Clone
        inline autoPtr<CompactListList<T, Container>> clone() const;


    // Member Functions

        // Access

            //- Return the primary size, i.e. the number of rows
            inline label size() const;

            //- Return true if the number of rows is zero
            inline bool empty() const;

            //- Return the offset table (= size()+1)
            inline const List<label>& offsets() const;

            //- Return non-const access to the offset table
            inline List<label>& offsets();

            //- Return the packed matrix of data
            inline const List<T>& m() const;

            //- Return non-const access to the packed matrix of data
            inline List<T>& m();


        // Edit

            //- Reset size of CompactListList.
            //  This form only allows contraction of the CompactListList.
            void setSize(const label mRows);

            //- Reset size of CompactListList.
            void setSize(const label mRows, const label nData);

            //- Reset sizes of CompactListList and value for new elements.
            void setSize(const label mRows, const label nData, const T&);

            //- Reset size of CompactListList.
            void setSize(const labelUList& rowSizes);

            //- Reset size of CompactListList.
            //  This form only allows contraction of the CompactListList.
            inline void resize(const label mRows);

            //- Reset size of CompactListList.
            inline void resize(const label mRows, const label nData);

            //- Reset sizes of CompactListList and value for new elements.
            inline void resize(const label mRows, const label nData, const T&);

            //- Reset size of CompactListList.
            inline void resize(const labelUList& rowSizes);

            //- Clear the CompactListList, i.e. set sizes to zero.
            void clear();

            //- Return sizes (to be used e.g. for construction)
            labelList sizes() const;

            //- Transfer the contents of the argument CompactListList
            //  into this CompactListList and annul the argument list.
            void transfer(CompactListList<T, Container>&);


        // Other

            //- Return index into m
            inline label index(const label row, const label col) const;

            //- Get row for index into m.
            inline label whichRow(const label index) const;

            //- Get column index (j) given above row
            inline label whichColumn(const label row, const label index) const;


    // Member Operators

        //- Return subscript-checked row as UList.
        inline UList<T> operator[](const label i);

        //- Return const subscript-checked row as UList.
        inline const UList<T> operator[](const label i) const;

        //- Return subscript-checked element.
        inline T& operator()(const label i, const label j);

        //- Return const subscript-checked element.
        inline const T& operator()(const label i, const label j) const;

        //- Return as List<Container>
        List<Container> operator()() const;

        //- Move assignment operator
        inline void operator=(CompactListList<T, Container>&&);

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // Istream operator

        //- Read CompactListList from Istream, discarding contents
        //  of existing CompactListList.
        friend Istream& operator>> <T, Container>
        (
            Istream&,
            CompactListList<T, Container>&
        );

        // Write CompactListList to Ostream.
        friend Ostream& operator<< <T, Container>
        (
            Ostream&,
            const CompactListList<T, Container>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "CompactListListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactListList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
