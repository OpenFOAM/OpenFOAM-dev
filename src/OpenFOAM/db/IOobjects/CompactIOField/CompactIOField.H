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
    Foam::CompactIOField

Description
    A Field of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOField except when binary output in
    case it writes a CompactListList.

    Useful for fields of small subfields e.g. in lagrangian

SourceFiles
    CompactIOField.C

\*---------------------------------------------------------------------------*/

#ifndef CompactIOField_H
#define CompactIOField_H

#include "IOField.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type, class BaseType> class CompactIOField;

template<class Type, class BaseType> Istream& operator>>
(
    Istream&,
    CompactIOField<Type, BaseType>&
);

template<class Type, class BaseType> Ostream& operator<<
(
    Ostream&,
    const CompactIOField<Type, BaseType>&
);

/*---------------------------------------------------------------------------*\
                       Class CompactIOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class BaseType>
class CompactIOField
:
    public regIOobject,
    public Field<Type>
{
    // Private Member Functions

        //- Read according to header type
        void readFromStream(const bool read = true);

public:

    //- Runtime type information
    TypeName("FieldField");


    // Constructors

        //- Construct from IOobject
        CompactIOField(const IOobject&);

        //- Construct from IOobject; does local processor require reading?
        CompactIOField(const IOobject&, const bool read);

        //- Construct from IOobject and size
        CompactIOField(const IOobject&, const label);

        //- Construct from IOobject and a Field
        CompactIOField(const IOobject&, const Field<Type>&);

        //- Move construct by transferring the Field contents
        CompactIOField(const IOobject&, Field<Type>&&);

        //- Move constructor
        CompactIOField(CompactIOField<Type, BaseType>&&);


    // Destructor

        virtual ~CompactIOField();


    // Member Functions

        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member Operators

        void operator=(const CompactIOField<Type, BaseType>&);
        void operator=(CompactIOField<Type, BaseType>&&);

        void operator=(const Field<Type>&);
        void operator=(Field<Type>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
