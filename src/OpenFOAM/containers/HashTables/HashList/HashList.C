/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "HashList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type, class Key, class Hash>
    const Key HashList<Type, Key, Hash>::nullKey = Key::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class Key, class Hash>
Foam::HashList<Type, Key, Hash>::HashList(const label size)
:
    List<Tuple2<Key, Type>>(size, Tuple2<Key, Type>(nullKey, Type()))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, class Key, class Hash>
inline Foam::label Foam::HashList<Type, Key, Hash>::capacity() const
{
    return List<Tuple2<Key, Type>>::size();
}


template<class Type, class Key, class Hash>
void Foam::HashList<Type, Key, Hash>::clear()
{
    List<Tuple2<Key, Type>>::operator=
    (
        Tuple2<Key, Type>(nullKey, Type())
    );
}


template<class Type, class Key, class Hash>
void Foam::HashList<Type, Key, Hash>::resizeAndClear(const label newSize)
{
    List<Tuple2<Key, Type>>::resize
    (
        newSize,
        Tuple2<Key, Type>(nullKey, Type())
    );
}


template<class Type, class Key, class Hash>
bool Foam::HashList<Type, Key, Hash>::insert(const Key& k, const Type& t)
{
    List<Tuple2<Key, Type>>& map = *this;

    const label n = map.size();

    const unsigned h = Hash()(k);

    for (label i = 0; i < n; i ++)
    {
        const label hi = (h + i) % n;

        if (map[hi].first() == nullKey)
        {
            map[hi] = Tuple2<Key, Type>(k, t);
            return true;
        }

        if (map[hi].first() == k)
        {
            return false;
        }
    }

    FatalErrorInFunction
        << "Hash list is full"
        << exit(FatalError);

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type, class Key, class Hash>
const Type& Foam::HashList<Type, Key, Hash>::operator[](const Key& k) const
{
    const List<Tuple2<Key, Type>>& map = *this;

    const label n = map.size();

    const unsigned h = Hash()(k);

    for (label i = 0; i < n; i ++)
    {
        const label hi = (h + i) % n;

        if (map[hi].first() == k)
        {
            return map[hi].second();
        }
    }

    FatalErrorInFunction
        << "Hash list does not contain key \"" << k << "\""
        << exit(FatalError);

    return NullObjectRef<Type>();
}


// ************************************************************************* //
