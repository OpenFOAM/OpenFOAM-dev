/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
template<class T, class BaseType> class CompactIOField;

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    CompactIOField<T, BaseType>&
);

template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const CompactIOField<T, BaseType>&
);

/*---------------------------------------------------------------------------*\
                        Class CompactIOField Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class CompactIOField
:
    public regIOobject,
    public Field<T>
{
    // Private Member Functions

        //- Read according to header type
        void readFromStream();

public:

    //- Runtime type information
    TypeName("FieldField");


    // Constructors

        //- Construct from IOobject
        CompactIOField(const IOobject&);

        //- Construct from IOobject and size
        CompactIOField(const IOobject&, const label);

        //- Construct from IOobject and a Field
        CompactIOField(const IOobject&, const Field<T>&);

        //- Construct by transferring the Field contents
        CompactIOField(const IOobject&, const Xfer<Field<T> >&);


    // Destructor

        virtual ~CompactIOField();


    // Member functions

        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member operators

        void operator=(const CompactIOField<T, BaseType>&);

        void operator=(const Field<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CompactIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
