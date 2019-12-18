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

\*---------------------------------------------------------------------------*/

#include "IOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const bool read)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        Istream& is = readStream(typeName, read);

        if (read)
        {
            is >> *this;
        }
        close();
    }
    else if (io.readOpt() == IOobject::READ_IF_PRESENT)
    {
        bool haveFile = headerOk();

        Istream& is = readStream(typeName, haveFile && read);

        if (read && haveFile)
        {
            is >> *this;
        }
        close();
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const label size)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        Field<Type>::setSize(size);
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const Field<Type>& f)
:
    regIOobject(io),
    Field<Type>(f)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, Field<Type>&& f)
:
    regIOobject(io),
    Field<Type>(move(f))
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const tmp<Field<Type>>& f)
:
    regIOobject(io),
    Field<Type>(f)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


template<class Type>
Foam::IOField<Type>::IOField(IOField<Type>&& f)
:
    regIOobject(move(f)),
    Field<Type>(move(f))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::IOField<Type>::~IOField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOField<Type>::writeData(Ostream& os) const
{
    return (os << static_cast<const Field<Type>&>(*this)).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::IOField<Type>::operator=(const IOField<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


template<class Type>
void Foam::IOField<Type>::operator=(IOField<Type>&& rhs)
{
    Field<Type>::operator=(move(rhs));
}


template<class Type>
void Foam::IOField<Type>::operator=(const Field<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


template<class Type>
void Foam::IOField<Type>::operator=(Field<Type>&& rhs)
{
    Field<Type>::operator=(move(rhs));
}


template<class Type>
void Foam::IOField<Type>::operator=(const tmp<Field<Type>>& rhs)
{
    Field<Type>::operator=(rhs);
}


// ************************************************************************* //
