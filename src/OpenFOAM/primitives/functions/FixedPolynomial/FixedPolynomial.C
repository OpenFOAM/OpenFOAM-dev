/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "FixedPolynomial.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN>::FixedPolynomial()
:
    VectorSpace<FixedPolynomial<T, PN>, T, PN>()
{}


template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN>::FixedPolynomial(const zero)
:
    VectorSpace<FixedPolynomial<T, PN>, T, PN>(Zero)
{}


template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN>::FixedPolynomial(const T (&coeffs)[PN])
:
    VectorSpace<FixedPolynomial<T, PN>, T, PN>()
{
    for (uint8_t i = 0; i < PN; ++ i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN>::FixedPolynomial(const UList<T>& coeffs)
:
    VectorSpace<FixedPolynomial<T, PN>, T, PN>()
{
    if (coeffs.size() != PN)
    {
        FatalErrorInFunction
            << "Given " << coeffs.size() << " coefficients for polynomial with "
            << PN << " terms" << exit(FatalError);
    }

    for (uint8_t i = 0; i < PN; ++ i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN>::FixedPolynomial(Istream& is)
:
    VectorSpace<FixedPolynomial<T, PN>, T, PN>(is)
{}


template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN>::FixedPolynomial
(
    const char* symbols,
    const Function1s::unitSets& units,
    const dictionary& dict
)
:
    VectorSpace<FixedPolynomial<T, PN>, T, PN>()
{
    unitSet unitsXpow(dimless);
    for (uint8_t i = 0; i < PN; ++ i)
    {
        this->v_[i] = dict.lookup<T>(word(symbols[i]), units.value/unitsXpow);
        unitsXpow.reset(units.x*unitsXpow);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, uint8_t PN>
T Foam::FixedPolynomial<T, PN>::value(const scalar x) const
{
    T result = this->v_[PN - 1];
    for (uint8_t i = 1; i < PN; ++ i)
    {
        result *= x;
        result += this->v_[PN - 1 - i];
    }
    return result;
}


template<class T, uint8_t PN>
T Foam::FixedPolynomial<T, PN>::derivative(const scalar x) const
{
    T result = this->v_[PN - 1]*(PN - 1);
    for (uint8_t i = 1; i < PN - 1; ++ i)
    {
        result *= x;
        result += this->v_[PN - 1 - i]*(PN - 1 - i);
    }
    return result;
}


template<class T, uint8_t PN>
T Foam::FixedPolynomial<T, PN>::integral(const scalar x) const
{
    T result = pTraits<T>::zero;
    for (uint8_t i = 0; i < PN; ++ i)
    {
        result += this->v_[PN - 1 - i]/(PN - i);
        result *= x;
    }
    return result;
}


template<class T, uint8_t PN>
T Foam::FixedPolynomial<T, PN>::integral
(
    const scalar x0,
    const scalar x1
) const
{
    return integral(x1) - integral(x0);
}


template<class T, uint8_t PN>
Foam::FixedPolynomial<T, PN + 1>
Foam::FixedPolynomial<T, PN>::integralPolynomial(const scalar x) const
{
    FixedPolynomial<T, PN + 1> result(Zero);
    for (uint8_t i = 0; i < PN; ++ i)
    {
        result[i + 1] = this->v_[i]/(i + 1);
    }
    result[0] = - result.value(x);
    return result;
}


template<class T, uint8_t PN>
void Foam::FixedPolynomial<T, PN>::write
(
    const char* symbols,
    const Function1s::unitSets& units,
    Ostream& os
) const
{
    FixedPolynomial<T, PN> copy(*this);
    convert(copy, units, unitSet::makeUserOp());
    for (uint8_t i = 0; i < PN; ++ i)
    {
        writeEntry(os, word(symbols[i]), copy[i]);
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T, uint8_t PN>
void Foam::convert
(
    FixedPolynomial<T, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeStandardOp&
)
{
    unitSet unitsXpow(dimless);
    for (uint8_t i = 0; i < PN; ++ i)
    {
        poly[i] = units.value.toStandard(unitsXpow.toUser(poly[i]));
        unitsXpow.reset(units.x*unitsXpow);
    }
}


template<class T, uint8_t PN>
void Foam::convert
(
    FixedPolynomial<T, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeUserOp&
)
{
    unitSet unitsXpow(dimless);
    for (uint8_t i = 0; i < PN; ++ i)
    {
        poly[i] = units.value.toUser(unitsXpow.toStandard(poly[i]));
        unitsXpow.reset(units.x*unitsXpow);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, uint8_t PN>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<FixedPolynomial<T, PN>>& ip
)
{
    const FixedPolynomial<T, PN>& poly = ip.t_;

    if (PN > 0)
    {
        os << poly[0];
    }
    for (label p = 1; p < PN; ++ p)
    {
        os << " + " << poly[p] << "*x";
        if (p > 1) os << "^" << p;
    }

    return os;
}


// ************************************************************************* //
