/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "GlobalIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::GlobalIOListBase<Container, IOContainer, Type>::GlobalIOListBase
(
    const IOobject& io
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOContainer<Type>>();

    readHeaderOk(IOstream::BINARY, IOContainer<Type>::typeName);
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::GlobalIOListBase<Container, IOContainer, Type>::GlobalIOListBase
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOContainer<Type>>();

    if (!readHeaderOk(IOstream::BINARY, IOContainer<Type>::typeName))
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
Foam::GlobalIOListBase<Container, IOContainer, Type>::GlobalIOListBase
(
    const IOobject& io,
    const Container<Type>& f
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOContainer<Type>>();

    if (!readHeaderOk(IOstream::BINARY, IOContainer<Type>::typeName))
    {
        Container<Type>::operator=(f);
    }
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::GlobalIOListBase<Container, IOContainer, Type>::GlobalIOListBase
(
    const IOobject& io,
    Container<Type>&& f
)
:
    regIOobject(io),
    Container<Type>(move(f))

{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOContainer<Type>>();

    readHeaderOk(IOstream::BINARY, IOContainer<Type>::typeName);
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::GlobalIOListBase<Container, IOContainer, Type>::GlobalIOListBase
(
    const GlobalIOListBase<Container, IOContainer, Type>& field
)
:
    regIOobject(field),
    Container<Type>(field)
{}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::GlobalIOListBase<Container, IOContainer, Type>::GlobalIOListBase
(
    GlobalIOListBase<Container, IOContainer, Type>&& field
)
:
    regIOobject(move(field)),
    Container<Type>(move(field))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
Foam::GlobalIOListBase<Container, IOContainer, Type>::~GlobalIOListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
bool Foam::GlobalIOListBase<Container, IOContainer, Type>::readData
(
    Istream& is
)
{
    is >> *this;
    return is.good();
}


template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
bool Foam::GlobalIOListBase<Container, IOContainer, Type>::writeData
(
    Ostream& os
) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
void Foam::GlobalIOListBase<Container, IOContainer, Type>::operator=
(
    const GlobalIOListBase<Container, IOContainer, Type>& rhs
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
void Foam::GlobalIOListBase<Container, IOContainer, Type>::operator=
(
    GlobalIOListBase<Container, IOContainer, Type>&& rhs
)
{
    Container<Type>::operator=(move(rhs));
}


// ************************************************************************* //
