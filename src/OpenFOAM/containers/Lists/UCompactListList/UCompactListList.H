/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::UCompactListList

Description
    Unallocated base class of CompactListList

SourceFiles
    UCompactListList.C
    UCompactListListI.H
    UCompactListListIO.C

\*---------------------------------------------------------------------------*/

#ifndef UCompactListList_H
#define UCompactListList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class UCompactListList;

template<class T> Ostream& operator<<(Ostream&, const UCompactListList<T>&);

/*---------------------------------------------------------------------------*\
                       Class UCompactListList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UCompactListList
{
    // Private Data

        //- Offset table
        UList<label> offsets_;

        //- Packed matrix of data
        UList<T> m_;


public:

    // Static Member Functions

        //- Return a null UCompactListList
        inline static const UCompactListList<T>& null();

        //- Set the sizes to match that of the given list-list
        template<class ListType, class ListListType>
        static void setSizeToListList
        (
            ListType& offsets,
            ListType& m,
            const ListListType& ll
        );

        //- Set the sizes and values to match that of the given list-list
        template<class ListType, class ListListType>
        static void setSizeAndValuesToListList
        (
            ListType& offsets,
            ListType& m,
            const ListListType& ll
        );


    // Constructors

        //- Null constructor.
        inline UCompactListList();

        //- Construct from components
        inline UCompactListList(const UList<label>& offsets, const UList<T>& m);


    // Member Functions

        // Access

            //- Return the primary size, i.e. the number of rows
            inline label size() const;

            //- Return true if the number of rows is zero
            inline bool empty() const;

            //- Return the offset table (= size()+1)
            inline const UList<label>& offsets() const;

            //- Return non-const access to the offset table
            inline UList<label>& offsets();

            //- Return the packed matrix of data
            inline const UList<T>& m() const;

            //- Return non-const access to the packed matrix of data
            inline UList<T>& m();


        //- Copy the ULists, but not the underlying data
        inline void shallowCopy(const UCompactListList<T>&);

        //- Copy the underlying data
        inline void deepCopy(const UCompactListList<T>&);

        //- Return index into m
        inline label index(const label row, const label col) const;

        //- Get row for index into m.
        inline label whichRow(const label index) const;

        //- Get column index (j) given above row
        inline label whichColumn(const label row, const label index) const;

        //- Return sizes (to be used e.g. for construction)
        labelList sizes() const;

        //- Convert to List<Container>
        template<class Container = List<T>>
        List<Container> list() const;


    // Member Operators

        //- Return subscript-checked row as UList.
        inline UList<T> operator[](const label i);

        //- Return const subscript-checked row as UList.
        inline const UList<T> operator[](const label i) const;

        //- Return subscript-checked element.
        inline T& operator()(const label i, const label j);

        //- Return const subscript-checked element.
        inline const T& operator()(const label i, const label j) const;

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // Istream operator

        // Write UCompactListList to Ostream.
        friend Ostream& operator<< <T>
        (
            Ostream&,
            const UCompactListList<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UCompactListListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UCompactListList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
