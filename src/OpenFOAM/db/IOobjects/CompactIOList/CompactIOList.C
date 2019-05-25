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

#include "CompactIOList.H"
#include "IOList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::readFromStream()
{
    Istream& is = readStream(word::null);

    if (headerClassName() == IOList<T>::typeName)
    {
        is >> static_cast<ListCompactIO<T, BaseType>&>(*this);
        close();
    }
    else if (headerClassName() == typeName)
    {
        is >> *this;
        close();
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "unexpected class name " << headerClassName()
            << " expected " << typeName << " or " << IOList<T>::typeName
            << endl
            << "    while reading object " << name()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList(const IOobject& io)
:
    regIOobject(io)
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromStream();
    }
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io)
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromStream();
    }
    else
    {
        this->setSize(size);
    }
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    const List<T>& list
)
:
    regIOobject(io)
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromStream();
    }
    else
    {
        ListCompactIO<T, BaseType>::operator=(list);
    }
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    List<T>&& list
)
:
    regIOobject(io),
    ListCompactIO<T, BaseType>(move(list))
{
    if
    (
        io.readOpt() == IOobject::MUST_READ
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromStream();
    }
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    CompactIOList<T, BaseType>&& list
)
:
    regIOobject(move(list)),
    List<T>(move(list))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::~CompactIOList()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    if (fmt == IOstream::ASCII)
    {
        // Change type to be non-compact format type
        const word oldTypeName = typeName;

        const_cast<word&>(typeName) = IOList<T>::typeName;

        bool good = regIOobject::writeObject(fmt, ver, cmp, write);

        // Change type back
        const_cast<word&>(typeName) = oldTypeName;

        return good;
    }
    else if (this->overflows())
    {
        WarningInFunction
            << "Overall number of elements of CompactIOList of size "
            << this->size() << " overflows the representation of a label"
            << endl << "    Switching to ascii writing" << endl;

        // Change type to be non-compact format type
        const word oldTypeName = typeName;

        const_cast<word&>(typeName) = IOList<T>::typeName;

        bool good = regIOobject::writeObject(IOstream::ASCII, ver, cmp, write);

        // Change type back
        const_cast<word&>(typeName) = oldTypeName;

        return good;
    }
    else
    {
        return regIOobject::writeObject(fmt, ver, cmp, write);
    }
}


template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=
(
    const CompactIOList<T, BaseType>& rhs
)
{
    ListCompactIO<T, BaseType>::operator=(rhs);
}


template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=
(
    CompactIOList<T, BaseType>&& rhs
)
{
    ListCompactIO<T, BaseType>::operator=(move(rhs));
}


template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=
(
    const List<T>& rhs
)
{
    ListCompactIO<T, BaseType>::operator=(rhs);
}


template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=
(
    List<T>&& rhs
)
{
    ListCompactIO<T, BaseType>::operator=(move(rhs));
}


// ************************************************************************* //
