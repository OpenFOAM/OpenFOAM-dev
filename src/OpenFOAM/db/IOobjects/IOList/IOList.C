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

#include "IOList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::IOList<T>::IOList(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOList::IOList(const IOobject&)")
            << "IOList " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOList does not support automatic rereading."
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


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const label size)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOList::IOList(const IOobject&, const label)")
            << "IOList " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOList does not support automatic rereading."
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
        List<T>::setSize(size);
    }
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const List<T>& list)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOList::IOList(const IOobject&, const List<T>&)")
            << "IOList " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOList does not support automatic rereading."
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
        List<T>::operator=(list);
    }
}


template<class T>
Foam::IOList<T>::IOList(const IOobject& io, const Xfer<List<T> >& list)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
            "IOList::IOList(const IOobject&, const Xfer<List<T> >&)"
        )   << "IOList " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOList does not support automatic rereading."
            << endl;
    }

    List<T>::transfer(list());

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

template<class T>
Foam::IOList<T>::~IOList()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOList<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOList<T>::operator=(const IOList<T>& rhs)
{
    List<T>::operator=(rhs);
}


template<class T>
void Foam::IOList<T>::operator=(const List<T>& rhs)
{
    List<T>::operator=(rhs);
}


// ************************************************************************* //
