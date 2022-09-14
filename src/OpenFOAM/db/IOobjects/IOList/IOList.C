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

#include "IOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    if (!IOContainer<Type>::rereading)
    {
        warnNoRereading<IOContainer<Type>>();
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
        // For if MUST_READ_IF_MODIFIED
        if (IOContainer<Type>::rereading)
        {
            addWatch();
        }

        readStream(IOContainer<Type>::typeName) >> *this;
        close();
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase
(
    const IOobject& io,
    const bool read
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    if (!IOContainer<Type>::rereading)
    {
        warnNoRereading<IOContainer<Type>>();
    }

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        // For if MUST_READ_IF_MODIFIED
        if (IOContainer<Type>::rereading)
        {
            addWatch();
        }

        Istream& is = readStream(IOContainer<Type>::typeName, read);

        if (read)
        {
            is >> *this;
        }

        close();
    }
    else if (io.readOpt() == IOobject::READ_IF_PRESENT)
    {
        bool haveFile = headerOk();

        Istream& is = readStream(IOContainer<Type>::typeName, haveFile && read);

        if (read && haveFile)
        {
            is >> *this;
        }

        close();
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    if (!IOContainer<Type>::rereading)
    {
        warnNoRereading<IOContainer<Type>>();
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
        // For if MUST_READ_IF_MODIFIED
        if (IOContainer<Type>::rereading)
        {
            addWatch();
        }

        readStream(IOContainer<Type>::typeName) >> *this;
        close();
    }
    else
    {
        Container<Type>::setSize(size);
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase
(
    const IOobject& io,
    const Container<Type>& l
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    if (!IOContainer<Type>::rereading)
    {
        warnNoRereading<IOContainer<Type>>();
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
        // For if MUST_READ_IF_MODIFIED
        if (IOContainer<Type>::rereading)
        {
            addWatch();
        }

        readStream(IOContainer<Type>::typeName) >> *this;
        close();
    }
    else
    {
        Container<Type>::operator=(l);
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase
(
    const IOobject& io,
    Container<Type>&& l
)
:
    regIOobject(io),
    Container<Type>(move(l))
{
    // Check for MUST_READ_IF_MODIFIED
    if (!IOContainer<Type>::rereading)
    {
        warnNoRereading<IOContainer<Type>>();
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
        // For if MUST_READ_IF_MODIFIED
        if (IOContainer<Type>::rereading)
        {
            addWatch();
        }

        readStream(IOContainer<Type>::typeName) >> *this;
        close();
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase
(
    const IOListBase<Container, IOContainer, Type>& f
)
:
    regIOobject(f),
    Container<Type>(f)
{}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::IOListBase
(
    IOListBase<Container, IOContainer, Type>&& f
)
:
    regIOobject(move(f)),
    Container<Type>(move(f))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::IOListBase<Container, IOContainer, Type>::~IOListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
bool Foam::IOListBase<Container, IOContainer, Type>::writeData
(
    Ostream& os
) const
{
    return (os << static_cast<const Container<Type>&>(*this)).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
void Foam::IOListBase<Container, IOContainer, Type>::operator=
(
    const IOListBase<Container, IOContainer, Type>& rhs
)
{
    Container<Type>::operator=(rhs);
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
void Foam::IOListBase<Container, IOContainer, Type>::operator=
(
    IOListBase<Container, IOContainer, Type>&& rhs
)
{
    Container<Type>::operator=(move(rhs));
}


// ************************************************************************* //
