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
        inline SubField(const SubList<Type>&);

        //- Construct from a UList\<Type\>, using the entire size
        explicit inline SubField(const UList<Type>&);

        //- Construct from a UList\<Type\> with a given size
        inline SubField
        (
            const UList<Type>& list,
            const label subSize
        );

        //- Construct from a UList\<Type\> with a given size and start index
        inline SubField
        (
            const UList<Type>& list,
            const label subSize,
            const label startIndex
        );

        //- Construct as copy
        inline SubField(const SubField<Type>&);


    // Member functions

        //- Return a null SubField
        static inline const SubField<Type>& null();

        //- Return a component field of the field
        inline tmp<Field<cmptType>> component(const direction) const;

        //- Return the field transpose (only defined for second rank tensors)
        tmp<Field<Type>> T() const;


    // Member operators

        //- Assignment via UList operator. Takes linear time.
        inline void operator=(const SubField<Type>&);

        //- Assignment via UList operator. Takes linear time.
        inline void operator=(const Field<Type>&);

        //- Assignment via UList operator. Takes linear time.
        template<class Form, direction Ncmpts>
        inline void operator=(const VectorSpace<Form, Type, Ncmpts>&);

        //- Allow cast to a const Field\<Type\>&
        inline operator const Field<Type>&() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
