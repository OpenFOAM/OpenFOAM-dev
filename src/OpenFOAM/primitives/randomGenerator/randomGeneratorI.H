/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "randomGenerator.H"
#include "Hash.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline uint64_t Foam::randomGenerator::sample()
{
    x_ = (A*x_ + C) % M;
    return x_ >> 17;
}


// * * * * * * * * * * * * * Protected Constructors  * * * * * * * * * * * * //

inline Foam::randomGenerator::randomGenerator(const pTraits<uint64_t> s)
:
    x_((s << 16) + 0x330E)
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline uint64_t Foam::randomGenerator::x() const
{
    return x_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::randomGenerator::randomGenerator(const label s)
:
    randomGenerator(pTraits<uint64_t>(s))
{}


inline Foam::randomGenerator::randomGenerator(const word& s)
:
    randomGenerator(pTraits<uint64_t>(Hash<word>()(s)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::randomGenerator::~randomGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::randomGenerator::scalar01()
{
    return scalar(sample())/(M >> 17);
}


inline Foam::tmp<Foam::scalarField> Foam::randomGenerator::scalar01
(
    const label n
)
{
    tmp<scalarField> tValues(new scalarField(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = scalar01();
    }
    return tValues;
}


inline Foam::scalar Foam::randomGenerator::scalarAB
(
    const scalar a,
    const scalar b
)
{
    return a + scalar01()*(b - a);
}


inline Foam::tmp<Foam::scalarField> Foam::randomGenerator::scalarAB
(
    const label n,
    const scalar a,
    const scalar b
)
{
    tmp<scalarField> tValues(new scalarField(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = scalarAB(a, b);
    }
    return tValues;
}


template<class Type>
inline Type Foam::randomGenerator::sample01()
{
    Type value;
    for (direction i = 0; i < pTraits<Type>::nComponents; ++ i)
    {
        value.component(i) = scalar01();
    }
    return value;
}


template<>
inline Foam::scalar Foam::randomGenerator::sample01()
{
    return scalar01();
}


template<>
inline Foam::label Foam::randomGenerator::sample01()
{
    return sample() % 2;
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::randomGenerator::sample01
(
    const label n
)
{
    tmp<Field<Type>> tValues(new Field<Type>(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = sample01<Type>();
    }
    return tValues;
}


template<class Type>
inline Type Foam::randomGenerator::sampleAB(const Type& a, const Type& b)
{
    return a + cmptMultiply(sample01<Type>(), b - a);
}


template<>
inline Foam::scalar Foam::randomGenerator::sampleAB
(
    const scalar& a,
    const scalar& b
)
{
    return scalarAB(a, b);
}


template<>
inline Foam::label Foam::randomGenerator::sampleAB
(
    const label& a,
    const label& b
)
{
    return a + sample() % (b - a);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::randomGenerator::sampleAB
(
    const label n,
    const Type& a,
    const Type& b
)
{
    tmp<Field<Type>> tValues(new Field<Type>(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = sampleAB(a, b);
    }
    return tValues;
}


template<class Container>
inline void Foam::randomGenerator::permute(Container& l)
{
    for (label i = 0; i < l.size(); ++ i)
    {
        Swap(l[i], l[sampleAB<label>(i, l.size())]);
    }
}


// ************************************************************************* //
