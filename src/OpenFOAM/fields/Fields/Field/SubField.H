/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::SubField

Description
    SubField is a Field obtained as a section of another Field.

    Thus it is itself unallocated so that no storage is allocated or
    deallocated during its use.  To achieve this behaviour, SubField is
    derived from a SubList rather than a List.

SourceFiles
    SubFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef SubField_H
#define SubField_H

#include "SubList.H"
#include "Field.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Pre-declare SubField and related Field type
template<class Type> class Field;
template<class Type> class SubField;

/*---------------------------------------------------------------------------*\
                          Class SubField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SubField
:
    public tmp<SubField<Type>>::refCount,
    public SubList<Type>
{

public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Constructors

        //- Construct from a SubList
        SubField(const SubList<Type>&);

        //- Construct from a UList\<Type\>, using the entire size
        explicit SubField(const UList<Type>&);

        //- Construct from a UList\<Type\> with a given size
        SubField(const UList<Type>& list, const label subSize);

        //- Construct from a UList\<Type\> with a given size and start index
        SubField
        (
            const UList<Type>& list,
            const label subSize,
            const label startIndex
        );

        //- Copy constructor
        SubField(const SubField<Type>&);

        //- Copy constructor or reuse as specified
        SubField(SubField<Type>&, bool reuse);


    // Static Member Functions

        //- Return a null SubField
        static inline const SubField<Type>& null()
        {
            return NullObjectRef<SubField<Type>>();
        }


    // Member Functions

        //- Return a component field of the field
        tmp<Field<cmptType>> component(const direction) const;

        //- Return the field transpose (only defined for second rank tensors)
        tmp<Field<Type>> T() const;


    // Member Operators

        void operator=(const SubField<Type>&);
        void operator=(const UList<Type>&);
        void operator=(const tmp<Field<Type>>&);
        void operator=(const Type&);
        void operator=(const zero);

        template<class Form, direction Ncmpts>
        void operator=(const VectorSpace<Form, Type, Ncmpts>&);

        void operator+=(const UList<Type>&);
        void operator+=(const tmp<Field<Type>>&);

        void operator-=(const UList<Type>&);
        void operator-=(const tmp<Field<Type>>&);

        void operator*=(const UList<scalar>&);
        void operator*=(const tmp<Field<scalar>>&);

        void operator/=(const UList<scalar>&);
        void operator/=(const tmp<Field<scalar>>&);

        void operator+=(const Type&);
        void operator-=(const Type&);

        void operator*=(const scalar&);
        void operator/=(const scalar&);

        //- Allow cast to a const Field\<Type\>&
        inline operator const Field<Type>&() const
        {
            return *reinterpret_cast<const Field<Type>*>(this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SubField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
