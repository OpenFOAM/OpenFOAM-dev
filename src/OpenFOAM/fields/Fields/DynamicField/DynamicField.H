/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::DynamicField

Description
    Dynamically sized Field.

SourceFiles
    DynamicFieldI.H
    DynamicField.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicField_H
#define DynamicField_H

#include "Field.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
class DynamicField;

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Ostream& operator<<
(
    Ostream&,
    const DynamicField<T, SizeInc, SizeMult, SizeDiv>&
);

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Istream& operator>>
(
    Istream&,
    DynamicField<T, SizeInc, SizeMult, SizeDiv>&
);


/*---------------------------------------------------------------------------*\
                        Class DynamicField Declaration
\*---------------------------------------------------------------------------*/

template<class T, unsigned SizeInc=0, unsigned SizeMult=2, unsigned SizeDiv=1>
class DynamicField
:
    public Field<T>
{
    static_assert
    (
        (SizeInc || SizeMult) && SizeDiv,
        "Avoid invalid sizing parameters"
    );

    // Private Data

        //- The capacity (allocated size) of the underlying field.
        label capacity_;


public:

    // Static Member Functions

        //- Return a null field
        inline static const DynamicField<T, SizeInc, SizeMult, SizeDiv>& null()
        {
            return *reinterpret_cast
            <
                DynamicField<T, SizeInc, SizeMult, SizeDiv>*
            >(0);
        }


    // Constructors

        //- Construct null
        inline DynamicField();

        //- Construct given size.
        explicit inline DynamicField(const label);

        //- Construct from UList. Size set to UList size.
        //  Also constructs from DynamicField with different sizing parameters.
        explicit inline DynamicField(const UList<T>&);

        //- Move constructor transferring the list contents
        explicit inline DynamicField(List<T>&&);

        //- Construct by 1 to 1 mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        inline DynamicField
        (
            const UList<T>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Copy constructor
        inline DynamicField(const DynamicField<T, SizeInc, SizeMult, SizeDiv>&);

        //- Move constructor
        inline DynamicField(DynamicField<T, SizeInc, SizeMult, SizeDiv>&&);

        //- Construct from Istream. Size set to size of list read.
        explicit DynamicField(Istream&);

        //- Clone
        tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv>> clone() const;


    // Member Functions

        // Access

            //- Size of the underlying storage.
            inline label capacity() const;

        // Edit

            //- Alter the size of the underlying storage.
            //  The addressed size will be truncated if needed to fit, but will
            //  remain otherwise untouched.
            //  Use this or reserve() in combination with append().
            inline void setCapacity(const label);

            //- Alter the addressed list size.
            //  New space will be allocated if required.
            //  Use this to resize the list prior to using the operator[] for
            //  setting values (as per List usage).
            inline void setSize(const label);

            //- Alter the addressed list size and fill new space with a
            //  constant.
            inline void setSize(const label, const T&);

            //- Alter the addressed list size.
            //  New space will be allocated if required.
            //  Use this to resize the list prior to using the operator[] for
            //  setting values (as per List usage).
            inline void resize(const label);

            //- Alter the addressed list size and fill new space with a
            //  constant.
            inline void resize(const label, const T&);

            //- Reserve allocation space for at least this size.
            //  Never shrinks the allocated size, use setCapacity() for that.
            inline void reserve(const label);

            //- Clear the addressed list, i.e. set the size to zero.
            //  Allocated size does not change
            inline void clear();

            //- Clear the list and delete storage.
            inline void clearStorage();

            //- Shrink the allocated space to the number of elements used.
            //  Returns a reference to the DynamicField.
            inline DynamicField<T, SizeInc, SizeMult, SizeDiv>& shrink();


        // Member Operators

            //- Append an element at the end of the list
            inline DynamicField<T, SizeInc, SizeMult, SizeDiv>& append
            (
                const T&
            );

            //- Append a List at the end of this list
            inline DynamicField<T, SizeInc, SizeMult, SizeDiv>& append
            (
                const UList<T>&
            );

            //- Remove and return the top element
            inline T remove();

            //- Return non-const access to an element, resizing list if
            //  necessary
            inline T& operator()(const label);

            //- Assignment of all addressed entries to the given value
            inline void operator=(const T&);

            //- Assignment operator
            inline void operator=
            (
                const DynamicField<T, SizeInc, SizeMult, SizeDiv>&
            );

            //- Move assignment operator
            inline void operator=
            (
                DynamicField<T, SizeInc, SizeMult, SizeDiv>&&
            );

            //- Assignment to UList
            inline void operator=(const UList<T>&);

            //- Move assignment to List
            inline void operator=(List<T>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DynamicField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
