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

\*---------------------------------------------------------------------------*/

#include "IOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOField::IOField(const IOobject&)")
            << "IOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOField does not support automatic rereading."
            << endl;
    }

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
Foam::IOField<Type>::IOField(const IOobject& io, const label size)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOField::IOField(const IOobject&, const label)")
            << "IOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOField does not support automatic rereading."
            << endl;
    }

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
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOField::IOField(const IOobject&, const Field<Type>&)")
            << "IOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOField does not support automatic rereading."
            << endl;
    }

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
        Field<Type>::operator=(f);
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const Xfer<Field<Type> >& f)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "IOField::IOField(const IOobject&, const Xfer<Field<Type> >&)"
        )   << "IOField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOField does not support automatic rereading."
            << endl;
    }

    Field<Type>::transfer(f());

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
void Foam::IOField<Type>::operator=(const Field<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


// ************************************************************************* //
