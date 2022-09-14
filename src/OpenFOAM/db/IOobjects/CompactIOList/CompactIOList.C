/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
bool
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
overflows() const
{
    label size = 0;

    forAll(*this, i)
    {
        label oldSize = size;

        size += this->operator[](i).size();

        if (size < oldSize)
        {
            return true;
        }
    }

    return false;
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
convertToCompact
(
    labelList& start,
    List<typename Type::value_type>& elems
) const
{
    start.setSize(this->size() + 1);

    start[0] = 0;

    for (label i = 1; i < start.size(); i++)
    {
        const label prev = start[i-1];

        start[i] = prev + this->operator[](i-1).size();

        if (start[i] < prev)
        {
            FatalErrorInFunction
                << "Overall number of elements " << start[i]
                << " of CompactIOListBase of size "
                << this->size() << " overflows the representation of a label"
                << endl << "Please recompile with a larger representation"
                << " for label" << exit(FatalError);
        }
    }

    elems.setSize(start[start.size() - 1]);

    label elemi = 0;

    forAll(*this, i)
    {
        const Type& subList = this->operator[](i);

        forAll(subList, j)
        {
            elems[elemi++] = subList[j];
        }
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
convertFromCompact
(
    const labelList& start,
    const List<typename Type::value_type>& elems
)
{
    this->setSize(start.size() - 1);

    forAll(*this, i)
    {
        Type& subList = this->operator[](i);

        label index = start[i];

        subList.setSize(start[i+1] - index);

        forAll(subList, j)
        {
            subList[j] = elems[index++];
        }
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
readFromStream(const bool read)
{
    Istream& is = readStream(word::null, read);

    if (read)
    {
        if (headerClassName() == IOContainer<Type>::typeName)
        {
            is >> static_cast<Container<Type>&>(*this);
            close();
        }
        else if (headerClassName() == CompactIOContainer<Type>::typeName)
        {
            is >> *this;
            close();
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "unexpected class name " << headerClassName()
                << " expected " << CompactIOContainer<Type>::typeName << " or "
                << IOContainer<Type>::typeName << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
CompactIOListBase(const IOobject& io)
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


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
CompactIOListBase
(
    const IOobject& io,
    const bool read
)
:
    regIOobject(io)
{
    if (io.readOpt() == IOobject::MUST_READ)
    {
        readFromStream(read);
    }
    else if (io.readOpt() == IOobject::READ_IF_PRESENT)
    {
        bool haveFile = headerOk();
        readFromStream(read && haveFile);
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
CompactIOListBase
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


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
CompactIOListBase
(
    const IOobject& io,
    const Container<Type>& l
)
:
    regIOobject(io),
    Container<Type>(l)
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


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
CompactIOListBase
(
    const IOobject& io,
    Container<Type>&& l
)
:
    regIOobject(io),
    Container<Type>(move(l))
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


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
CompactIOListBase
(
    CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>&& l
)
:
    regIOobject(move(l)),
    Container<Type>(move(l))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
~CompactIOListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
bool
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
writeObject
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
        const word oldTypeName = CompactIOContainer<Type>::typeName;

        const_cast<word&>(CompactIOContainer<Type>::typeName) =
            IOContainer<Type>::typeName;

        bool good = regIOobject::writeObject(fmt, ver, cmp, write);

        // Change type back
        const_cast<word&>(CompactIOContainer<Type>::typeName) = oldTypeName;

        return good;
    }
    else if (this->overflows())
    {
        WarningInFunction
            << "Overall number of elements of CompactIOListBase of size "
            << this->size() << " overflows the representation of a label"
            << endl << "    Switching to ascii writing" << endl;

        // Change type to be non-compact format type
        const word oldTypeName = CompactIOContainer<Type>::typeName;

        const_cast<word&>(CompactIOContainer<Type>::typeName) =
            IOContainer<Type>::typeName;

        bool good = regIOobject::writeObject(IOstream::ASCII, ver, cmp, write);

        // Change type back
        const_cast<word&>(CompactIOContainer<Type>::typeName) = oldTypeName;

        return good;
    }
    else
    {
        return regIOobject::writeObject(fmt, ver, cmp, write);
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
bool
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
operator=
(
    const CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>&
        rhs
)
{
    Container<Type>::operator=(rhs);
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void
Foam::CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>::
operator=
(
    CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>&& rhs
)
{
    Container<Type>::operator=(move(rhs));
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void Foam::writeEntry
(
    Ostream& os,
    const CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>& l
)
{
    os << l;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>& l
)
{
    // ASCII is not compacted
    if (is.format() == IOstream::ASCII)
    {
        is >> static_cast<Container<Type>&>(l);
    }
    else
    {
        labelList start(is);
        List<typename Type::value_type> elems(is);
        l.convertFromCompact(start, elems);
    }

    return is;
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>& l
)
{
    // ASCII is not compacted
    if (os.format() == IOstream::ASCII)
    {
        os << static_cast<const Container<Type>&>(l);
    }
    else
    {
        labelList start;
        List<typename Type::value_type> elems;
        l.convertToCompact(start, elems);
        os << start << elems;
    }

    return os;
}


// ************************************************************************* //
