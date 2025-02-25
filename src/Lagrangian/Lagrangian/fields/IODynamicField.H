/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::IODynamicField

Description
    A dynamic field of type \<Type\> with automated input and output.

\*---------------------------------------------------------------------------*/

#ifndef IODynamicField_H
#define IODynamicField_H

#include "IOList.H"
#include "DynamicField.H"
#include "OldTimeField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
using IODynamicFieldContainer = DynamicField<Type>;

template<class Type>
class IODynamicField;

template<class Type>
struct OldTimeField0Type<IODynamicField<Type>>
{
    typedef IODynamicField<Type> type;
};

template<class Type>
struct OldTimeOtherFieldType<IODynamicField<Type>>
{
    template<template<class> class OtherPrimitiveField>
    using type = IODynamicField<Type>;
};

/*---------------------------------------------------------------------------*\
                          Class IODynamicField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IODynamicField
:
    public IOListBase<IODynamicFieldContainer, IODynamicField, Type>,
    public OldTimeField<IODynamicField<Type>>
{
    // Private Type Definitions

        //- The base IO list class
        typedef IOListBase<IODynamicFieldContainer, IODynamicField, Type> Base;


public:

    //- Runtime type information
    TypeName("DynamicField");


    // Constructors

        //- Forward construction
        template<class ... Args>
        IODynamicField(Args ... args)
        :
            Base(args ...),
            OldTimeField<IODynamicField<Type>>(this->time().timeIndex())
        {}

        //- Copy construct
        IODynamicField(const IODynamicField& iodf)
        :
            Base(iodf),
            OldTimeField<IODynamicField<Type>>(iodf)
        {}

        //- Move construct
        IODynamicField(IODynamicField&& iodf)
        :
            Base(std::move(iodf)),
            OldTimeField<IODynamicField<Type>>(std::move(iodf))
        {}


    // Member Operators

        //- Copy-assign
        void operator=(const IODynamicField<Type>& iodf)
        {
            Base::operator=(iodf);
        }

        //- Move-assign
        void operator=(IODynamicField<Type>&& iodf)
        {
            Base::operator=(std::move(iodf));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
