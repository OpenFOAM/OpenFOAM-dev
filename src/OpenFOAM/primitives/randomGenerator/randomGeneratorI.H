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

inline uint64_t Foam::randomGenerator::seed::x(const bool global) const
{
    const uint64_t localS =
        global ? s_ : s_ + (UINT64_MAX/Pstream::nProcs())*Pstream::myProcNo();

    return (localS << 16) + 0x330E;
}


inline void Foam::randomGenerator::checkSync() const
{
    if (global_)
    {
        uint64_t xMaster = x_;
        Pstream::scatter(xMaster);
        if (xMaster != x_)
        {
            FatalErrorInFunction
                << "Global random number generator is not synchronised"
                << exit(FatalError);
        }
    }
}


inline uint64_t Foam::randomGenerator::sample()
{
    x_ = (A*x_ + C) % M;
    return x_ >> 17;
}


inline Foam::scalar Foam::randomGenerator::scalar01NoCheckSync()
{
    return scalar(sample())/(M >> 17);
}


inline Foam::scalar Foam::randomGenerator::scalarABNoCheckSync
(
    const scalar a,
    const scalar b
)
{
    return a + scalar01NoCheckSync()*(b - a);
}


template<class Type>
inline Type Foam::randomGenerator::sample01NoCheckSync()
{
    Type value;
    for (direction i = 0; i < pTraits<Type>::nComponents; ++ i)
    {
        value.component(i) = scalar01NoCheckSync();
    }
    return value;
}


template<>
inline Foam::scalar Foam::randomGenerator::sample01NoCheckSync()
{
    return scalar01NoCheckSync();
}


template<>
inline Foam::label Foam::randomGenerator::sample01NoCheckSync()
{
    return sample() % 2;
}


template<class Type>
inline Type Foam::randomGenerator::sampleABNoCheckSync
(
    const Type& a,
    const Type& b
)
{
    return a + cmptMultiply(sample01NoCheckSync<Type>(), b - a);
}


template<>
inline Foam::scalar Foam::randomGenerator::sampleABNoCheckSync
(
    const scalar& a,
    const scalar& b
)
{
    return scalarABNoCheckSync(a, b);
}


template<>
inline Foam::label Foam::randomGenerator::sampleABNoCheckSync
(
    const label& a,
    const label& b
)
{
    return a + sample() % (b - a);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::randomGenerator::seed::seed(const label s)
:
    s_(s)
{}


inline Foam::randomGenerator::seed::seed(const word& s)
:
    s_(Hash<word>()(s))
{}


inline Foam::randomGenerator::randomGenerator(const seed s, const bool global)
:
    global_(global),
    x_(s.x(global))
{
    checkSync();
}


inline Foam::randomGenerator::randomGenerator(const randomGenerator& rndGen)
:
    global_(rndGen.global_),
    x_(rndGen.x_)
{
    checkSync();
}


inline Foam::randomGenerator::randomGenerator(randomGenerator&& rndGen)
:
    randomGenerator(static_cast<const randomGenerator&>(rndGen))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::randomGenerator::~randomGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::randomGenerator::scalar01()
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    return scalar01NoCheckSync();
}


inline Foam::tmp<Foam::scalarField> Foam::randomGenerator::scalar01
(
    const label n
)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    tmp<scalarField> tValues(new scalarField(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = scalar01NoCheckSync();
    }
    return tValues;
}


inline Foam::scalar Foam::randomGenerator::scalarAB
(
    const scalar a,
    const scalar b
)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    return scalarABNoCheckSync(a, b);
}


inline Foam::tmp<Foam::scalarField> Foam::randomGenerator::scalarAB
(
    const label n,
    const scalar a,
    const scalar b
)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    tmp<scalarField> tValues(new scalarField(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = scalarABNoCheckSync(a, b);
    }
    return tValues;
}


template<class Type>
inline Type Foam::randomGenerator::sample01()
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    return sample01NoCheckSync<Type>();
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::randomGenerator::sample01
(
    const label n
)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    tmp<Field<Type>> tValues(new Field<Type>(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = sample01NoCheckSync<Type>();
    }
    return tValues;
}


template<class Type>
inline Type Foam::randomGenerator::sampleAB(const Type& a, const Type& b)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    return sampleABNoCheckSync(a, b);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::randomGenerator::sampleAB
(
    const label n,
    const Type& a,
    const Type& b
)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    tmp<Field<Type>> tValues(new Field<Type>(n));
    for (label i = 0; i < n; ++ i)
    {
        tValues.ref()[i] = sampleABNoCheckSync(a, b);
    }
    return tValues;
}


template<class Container>
inline void Foam::randomGenerator::permute(Container& l)
{
    #ifdef FULLDEBUG
    checkSync();
    #endif

    for (label i = 0; i < l.size(); ++ i)
    {
        Swap(l[i], l[sampleABNoCheckSync<label>(i, l.size())]);
    }
}


inline Foam::randomGenerator Foam::randomGenerator::generator()
{
    return randomGenerator(sample(), global_);
}


// ************************************************************************* //
