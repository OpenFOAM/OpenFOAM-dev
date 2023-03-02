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

SourceFiles
    CompactListList.C
    CompactListListI.H
    CompactListListIO.C

\*---------------------------------------------------------------------------*/

#ifndef CompactListList_H
#define CompactListList_H

#include "UCompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class CompactListList;

template<class T> Istream& operator>>(Istream&, CompactListList<T>&);


/*---------------------------------------------------------------------------*\
                       Class CompactListList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CompactListList
:
    public UCompactListList<T>
{
    // Private Data

        //- Offset table
        List<label> offsets_;

        //- Packed matrix of data
        List<T> m_;


public:

    // Static Member Functions

        //- Return a null CompactListList
        inline static const CompactListList<T>& null();


    // Constructors

        //- Null constructor.
        inline CompactListList();

        //- Construct from components
        inline CompactListList(const UList<label>& offsets, const UList<T>& m);

        //- Construct as copy or re-use as specified.
        inline CompactListList(CompactListList<T>&, bool reuse);

        //- Move constructor
        inline CompactListList(CompactListList<T>&&);

        //- Construct given size of offset table (number of rows),
        //  the number of data and a value for all elements.
        inline CompactListList(const label mRows, const label nData, const T&);

        //- Construct from UList<T2>
        template<class T2>
        CompactListList(const UList<T2>&);

        //- Construct from UIndirectList<T2>
        template<class T2>
        CompactListList(const UIndirectList<T2>&);

        //- Construct given list of row-sizes and a value for all elements
        CompactListList(const labelUList& rowSizes, const T&);

        //- Construct from Istream.
        CompactListList(Istream&);

        //- Clone
        inline autoPtr<CompactListList<T>> clone() const;


    // Member Functions

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

            //- Reset size of CompactListList to match the given UList<T2>
            template<class T2>
            void setSize(const UList<T2>&);

            //- Reset size of CompactListList to match the given
            //  UIndirectList<T2>
            template<class T2>
            void setSize(const UIndirectList<T2>&);

            //- Reset size of CompactListList.
            //  This form only allows contraction of the CompactListList.
            inline void resize(const label mRows);

            //- Reset size of CompactListList.
            inline void resize(const label mRows, const label nData);

            //- Reset sizes of CompactListList and value for new elements.
            inline void resize(const label mRows, const label nData, const T&);

            //- Reset size of CompactListList.
            inline void resize(const labelUList& rowSizes);

            //- Reset size of CompactListList to match the given
            //  UList<T2>
            template<class T2>
            void resize(const UList<T2>&);

            //- Reset size of CompactListList to match the given
            //  UIndirectList<T2>
            template<class T2>
            void resize(const UIndirectList<T2>&);

            //- Clear the CompactListList, i.e. set sizes to zero.
            void clear();

            //- Transfer the contents of the argument CompactListList
            //  into this CompactListList and annul the argument list.
            void transfer(CompactListList<T>&);


        //- Disallow shallowCopy
        void shallowCopy(const UCompactListList<T>&) = delete;


    // Istream operator

        //- Read CompactListList from Istream, discarding contents
        //  of existing CompactListList.
        friend Istream& operator>> <T>
        (
            Istream&,
            CompactListList<T>&
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
