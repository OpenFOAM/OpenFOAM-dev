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

#include "FixedLaurentPolynomial.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int8_t P0, int8_t PN>
Foam::FixedLaurentPolynomial<T, P0, PN>::FixedLaurentPolynomial()
:
    VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>()
{}


template<class T, int8_t P0, int8_t PN>
Foam::FixedLaurentPolynomial<T, P0, PN>::FixedLaurentPolynomial(const zero)
:
    VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>(Zero)
{}


template<class T, int8_t P0, int8_t PN>
Foam::FixedLaurentPolynomial<T, P0, PN>::FixedLaurentPolynomial
(
    const T (&coeffs)[PN - P0]
)
:
    VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>()
{
    for (int8_t i = 0; i < PN - P0; ++ i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<class T, int8_t P0, int8_t PN>
Foam::FixedLaurentPolynomial<T, P0, PN>::FixedLaurentPolynomial
(
    const UList<T>& coeffs
)
:
    VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>()
{
    if (coeffs.size() != PN - P0)
    {
        FatalErrorInFunction
            << "Given " << coeffs.size() << " coefficients for polynomial with "
            << PN - P0 << " terms" << exit(FatalError);
    }

    for (int8_t i = 0; i < PN - P0; ++ i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<class T, int8_t P0, int8_t PN>
Foam::FixedLaurentPolynomial<T, P0, PN>::FixedLaurentPolynomial(Istream& is)
:
    VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>(is)
{}


template<class T, int8_t P0, int8_t PN>
Foam::FixedLaurentPolynomial<T, P0, PN>::FixedLaurentPolynomial
(
    const char* symbols,
    const Function1s::unitSets& units,
    const dictionary& dict
)
:
    VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>()
{
    // Negative powers
    unitSet unitsXpowNeg(dimless);
    for (int8_t p = 0; p > PN; -- p)
    {
        unitsXpowNeg.reset(unitsXpowNeg/units.x);
    }
    for (int8_t p = min(PN, int8_t(0)) - 1; p > P0 - 1; -- p)
    {
        unitsXpowNeg.reset(unitsXpowNeg/units.x);
        this->v_[p - P0] =
            dict.lookup<T>(word(symbols[p - P0]), units.value/unitsXpowNeg);
    }

    // Positive powers and constant term
    unitSet unitsXpowPos(dimless);
    for (int8_t p = 0; p < P0; ++ p)
    {
        unitsXpowPos.reset(unitsXpowPos*units.x);
    }
    for (int8_t p = max(P0, int8_t(0)); p < PN; ++ p)
    {
        this->v_[p - P0] =
            dict.lookup<T>(word(symbols[p - P0]), units.value/unitsXpowPos);
        unitsXpowPos.reset(unitsXpowPos*units.x);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int8_t P0, int8_t PN>
T Foam::FixedLaurentPolynomial<T, P0, PN>::value(const scalar x) const
{
    // Negative powers
    T resultNeg = pTraits<T>::zero;
    for (int8_t p = P0; p < min(PN, int8_t(0)); ++ p)
    {
        resultNeg += this->v_[p - P0];
        resultNeg /= x;
    }
    for (int8_t p = PN; p < 0; ++ p)
    {
        resultNeg /= x;
    }

    // Positive powers and constant term
    T resultPos = PN > 0 ? this->v_[PN - 1 - P0] : pTraits<T>::zero;
    for (int8_t p = PN - 2; p >= max(P0, int8_t(0)); -- p)
    {
        resultPos *= x;
        resultPos += this->v_[p - P0];
    }
    for (int8_t p = P0 - 1; p >= 0; -- p)
    {
        resultPos *= x;
    }

    return resultNeg + resultPos;
}


template<class T, int8_t P0, int8_t PN>
T Foam::FixedLaurentPolynomial<T, P0, PN>::derivative(const scalar x) const
{
    // Negative powers
    T resultNeg = pTraits<T>::zero;
    for (int8_t p = P0; p <= min(PN, int8_t(0)); ++ p)
    {
        resultNeg += this->v_[p - P0]*p;
        resultNeg /= x;
    }
    for (int8_t p = PN; p < 0; ++ p)
    {
        resultNeg /= x;
    }

    // Positive powers
    T resultPos = PN > 0 ? this->v_[PN - 1 - P0]*(PN - 1) : pTraits<T>::zero;
    for (int8_t p = PN - 2; p >= max(P0, int8_t(1)); -- p)
    {
        resultPos *= x;
        resultPos += this->v_[p - P0]*p;
    }
    for (int8_t p = P0 - 2; p >= 0; -- p)
    {
        resultPos *= x;
    }

    return resultNeg + resultPos;
}


template<class T, int8_t P0, int8_t PN>
T Foam::FixedLaurentPolynomial<T, P0, PN>::integral(const scalar x) const
{
    // Negative powers
    T resultNeg = pTraits<T>::zero;
    for (int8_t p = P0; p < min(PN, int8_t(-1)); ++ p)
    {
        resultNeg /= x;
        resultNeg += this->v_[p - P0]/(p + 1);
    }
    for (int8_t p = PN; p < 0; ++ p)
    {
        resultNeg /= x;
    }

    // Log term
    if (P0 < 0 && 0 <= PN)
    {
        resultNeg /= x;
        resultNeg += this->v_[-1 - P0]*log(x);
    }

    // Positive powers
    T resultPos = pTraits<T>::zero;
    for (int8_t p = PN - 1; p >= max(P0, int8_t(0)); -- p)
    {
        resultPos += this->v_[p - P0]/(p + 1);
        resultPos *= x;
    }
    for (int8_t p = P0 - 1; p >= 0; -- p)
    {
        resultPos *= x;
    }

    return resultNeg + resultPos;
}


template<class T, int8_t P0, int8_t PN>
T Foam::FixedLaurentPolynomial<T, P0, PN>::integral
(
    const scalar x0,
    const scalar x1
) const
{
    return integral(x1) - integral(x0);
}


template<class T, int8_t P0, int8_t PN>
void Foam::FixedLaurentPolynomial<T, P0, PN>::write
(
    const char* symbols,
    const Function1s::unitSets& units,
    Ostream& os
) const
{
    FixedLaurentPolynomial<T, P0, PN> copy(*this);
    convert(copy, units, unitSet::makeUserOp());
    for (uint8_t i = 0; i < PN - P0; ++ i)
    {
        writeEntry(os, word(symbols[i]), copy[i]);
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T, int8_t P0, int8_t PN>
void Foam::convert
(
    FixedLaurentPolynomial<T, P0, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeStandardOp&
)
{
    // Negative powers
    unitSet unitsXpowNeg(dimless);
    for (int8_t p = 0; p > PN; -- p)
    {
        unitsXpowNeg.reset(unitsXpowNeg/units.x);
    }
    for (int8_t p = min(PN, int8_t(0)) - 1; p > P0 - 1; -- p)
    {
        unitsXpowNeg.reset(unitsXpowNeg/units.x);
        poly[p - P0] =
            units.value.toStandard(unitsXpowNeg.toUser(poly[p - P0]));
    }

    // Positive powers and constant term
    unitSet unitsXpowPos(dimless);
    for (int8_t p = 0; p < P0; ++ p)
    {
        unitsXpowPos.reset(unitsXpowPos*units.x);
    }
    for (int8_t p = max(P0, int8_t(0)); p < PN; ++ p)
    {
        poly[p - P0] =
            units.value.toStandard(unitsXpowPos.toUser(poly[p - P0]));
        unitsXpowPos.reset(unitsXpowPos*units.x);
    }
}


template<class T, int8_t P0, int8_t PN>
void Foam::convert
(
    FixedLaurentPolynomial<T, P0, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeUserOp&
)
{
    // Negative powers
    unitSet unitsXpowNeg(dimless);
    for (int8_t p = 0; p > PN; -- p)
    {
        unitsXpowNeg.reset(unitsXpowNeg/units.x);
    }
    for (int8_t p = min(PN, int8_t(0)) - 1; p > P0 - 1; -- p)
    {
        unitsXpowNeg.reset(unitsXpowNeg/units.x);
        poly[p - P0] =
            units.value.toUser(unitsXpowNeg.toStandard(poly[p - P0]));
    }

    // Positive powers and constant term
    unitSet unitsXpowPos(dimless);
    for (int8_t p = 0; p < P0; ++ p)
    {
        unitsXpowPos.reset(unitsXpowPos*units.x);
    }
    for (int8_t p = max(P0, int8_t(0)); p < PN; ++ p)
    {
        poly[p - P0] =
            units.value.toUser(unitsXpowPos.toStandard(poly[p - P0]));
        unitsXpowPos.reset(unitsXpowPos*units.x);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int8_t P0, int8_t PN>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<FixedLaurentPolynomial<T, P0, PN>>& ip
)
{
    const FixedLaurentPolynomial<T, P0, PN>& poly = ip.t_;

    List<string> terms(PN - P0);

    // Negative powers
    for (int8_t p = min(PN, int8_t(0)) - 1; p > P0 - 1; -- p)
    {
        terms[p - P0] = Foam::name(poly.v_[p - P0]) + "/x";
        if (p < -1) terms[p - P0] += "^" + Foam::name(- p);
    }

    // Positive powers and constant term
    for (int8_t p = max(P0, int8_t(0)); p < PN; ++ p)
    {
        terms[p - P0] = Foam::name(poly.v_[p - P0]);
        if (p > 0) terms[p - P0] += "*x";
        if (p > 1) terms[p - P0] += "^" + Foam::name(p);
    }

    // Print
    if (PN > P0)
    {
        os << terms[0].c_str();
    }
    for (int8_t i = 1; i < PN - P0; ++ i)
    {
        os << " + " << terms[i].c_str();
    }

    return os;
}


// ************************************************************************* //
