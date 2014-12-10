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

#include "CompactIOField.H"
#include "labelList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOField<T, BaseType>::readFromStream()
{
    Istream& is = readStream(word::null);

    if (headerClassName() == IOField<T>::typeName)
    {
        is >> static_cast<Field<T>&>(*this);
        close();
    }
    else if (headerClassName() == typeName)
    {
        is >> *this;
        close();
    }
    else
    {
        FatalIOErrorIn
        (
            "CompactIOField<T, BaseType>::readFromStream()",
            is
        )   << "unexpected class name " << headerClassName()
            << " expected " << typeName << " or " << IOField<T>::typeName
            << endl
            << "    while reading object " << name()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField(const IOobject& io)
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
Foam::CompactIOField<T, BaseType>::CompactIOField
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
        Field<T>::setSize(size);
    }
}


template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField
(
    const IOobject& io,
    const Field<T>& list
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
        Field<T>::operator=(list);
    }
}


template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField
(
    const IOobject& io,
    const Xfer<Field<T> >& list
)
:
    regIOobject(io)
{
    Field<T>::transfer(list());

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromStream();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::~CompactIOField()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOField<T, BaseType>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (fmt == IOstream::ASCII)
    {
        // Change type to be non-compact format type
        const word oldTypeName = typeName;

        const_cast<word&>(typeName) = IOField<T>::typeName;

        bool good = regIOobject::writeObject(fmt, ver, cmp);

        // Change type back
        const_cast<word&>(typeName) = oldTypeName;

        return good;
    }
    else
    {
        return regIOobject::writeObject(fmt, ver, cmp);
    }
}


template<class T, class BaseType>
bool Foam::CompactIOField<T, BaseType>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOField<T, BaseType>::operator=
(
    const CompactIOField<T, BaseType>& rhs
)
{
    Field<T>::operator=(rhs);
}


template<class T, class BaseType>
void Foam::CompactIOField<T, BaseType>::operator=(const Field<T>& rhs)
{
    Field<T>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::CompactIOField<T, BaseType>& L
)
{
    // Read compact
    const labelList start(is);
    const Field<BaseType> elems(is);

    // Convert
    L.setSize(start.size()-1);

    forAll(L, i)
    {
        T& subField = L[i];

        label index = start[i];
        subField.setSize(start[i+1] - index);

        forAll(subField, j)
        {
            subField[j] = elems[index++];
        }
    }

    return is;
}


template<class T, class BaseType>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::CompactIOField<T, BaseType>& L
)
{
    // Keep ascii writing same.
    if (os.format() == IOstream::ASCII)
    {
        os << static_cast<const Field<T>&>(L);
    }
    else
    {
        // Convert to compact format
        labelList start(L.size()+1);

        start[0] = 0;
        for (label i = 1; i < start.size(); i++)
        {
            start[i] = start[i-1]+L[i-1].size();
        }

        Field<BaseType> elems(start[start.size()-1]);

        label elemI = 0;
        forAll(L, i)
        {
            const T& subField = L[i];

            forAll(subField, j)
            {
                elems[elemI++] = subField[j];
            }
        }
        os << start << elems;
    }

    return os;
}


// ************************************************************************* //
