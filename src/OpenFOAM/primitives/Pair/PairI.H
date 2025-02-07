/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "Pair.H"
#include "tmp.H"
#include "token.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
inline const Foam::Pair<Type>& Foam::Pair<Type>::null()
{
    return NullObjectRef<Pair<Type>>();
}


template<class Type>
inline Foam::label Foam::Pair<Type>::size()
{
    return 2;
}


template<class Type>
inline void Foam::Pair<Type>::checkIndex(const label i)
{
    if (i < 0 || unsigned(i) > 1)
    {
        FatalErrorInFunction
            << "index " << i << " out of range 0 ... 1"
            << abort(FatalError);
    }
}


template<class Type>
inline int Foam::Pair<Type>::compare(const Pair<Type>& a, const Pair<Type>& b)
{
    if (a == b)
    {
        return 1;
    }
    else if (a == reverse(b))
    {
        return -1;
    }
    else
    {
        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::Pair<Type>::Pair()
{}


template<class Type>
inline Foam::Pair<Type>::Pair(const Type& f, const Type& s)
:
    f_(f),
    s_(s)
{}


template<class Type>
inline Foam::Pair<Type>::Pair(Istream& is)
{
    is >> *this;
}


template<class Type>
template<class HashT>
inline Foam::Pair<Type>::Hash<HashT>::Hash()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Type& Foam::Pair<Type>::first() const
{
    return f_;
}


template<class Type>
inline Type& Foam::Pair<Type>::first()
{
    return f_;
}


template<class Type>
inline const Type& Foam::Pair<Type>::second() const
{
    return s_;
}


template<class Type>
inline Type& Foam::Pair<Type>::second()
{
    return s_;
}


template<class Type>
inline const Type& Foam::Pair<Type>::other(const Type& a) const
{
    if (first() == second())
    {
        FatalErrorInFunction
            << "Call to other only valid for Pair with differing"
            << " elements:" << *this << abort(FatalError);
        return first();
    }
    else if (first() == a)
    {
        return second();
    }
    else
    {
        if (second() != a)
        {
            FatalErrorInFunction
                << "Pair " << *this
                << " does not contain " << a << abort(FatalError);
        }
        return first();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline Type& Foam::Pair<Type>::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return i ? s_ : f_;
}


template<class Type>
inline const Type& Foam::Pair<Type>::operator[](const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return i ? s_ : f_;
}


template<class Type>
template<class HashT>
inline unsigned Foam::Pair<Type>::Hash<HashT>::operator()
(
    const Pair<Type>& p,
    unsigned seed
) const
{
    return HashT()(p.second(), HashT()(p.first(), seed));
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::Pair<Type> Foam::reverse(const Pair<Type>& p)
{
    return Pair<Type>(p.second(), p.first());
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::operator==(const Pair<Type>& a, const Pair<Type>& b)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class Type>
inline bool Foam::operator!=(const Pair<Type>& a, const Pair<Type>& b)
{
    return !(a == b);
}


template<class Type>
inline bool Foam::operator<(const Pair<Type>& a, const Pair<Type>& b)
{
    return
    (
        a.first() < b.first()
     ||
        (
            !(b.first() < a.first())
         && a.second() < b.second()
        )
    );
}


template<class Type>
inline bool Foam::operator<=(const Pair<Type>& a, const Pair<Type>& b)
{
    return !(b < a);
}


template<class Type>
inline bool Foam::operator>(const Pair<Type>& a, const Pair<Type>& b)
{
    return (b < a);
}


template<class Type>
inline bool Foam::operator>=(const Pair<Type>& a, const Pair<Type>& b)
{
    return !(a < b);
}


template<class Type>
inline Foam::Istream& Foam::operator>>(Istream& is, Pair<Type>& p)
{
    is.readBegin("Pair");
    is >> p.first() >> p.second();
    is.readEnd("Pair");

    // Check state of Istream
    is.check("operator>>(Istream&, Pair<Type>&)");

    return is;
}


template<class Type>
inline Foam::Ostream& Foam::operator<<(Ostream& os, const Pair<Type>& p)
{
    os  << token::BEGIN_LIST
        << p.first() << token::SPACE << p.second()
        << token::END_LIST;

    return os;
}


template<class Type>
inline void Foam::writeEntry(Ostream& os, const Pair<Type>& p)
{
    os << p;
}


// ************************************************************************* //
