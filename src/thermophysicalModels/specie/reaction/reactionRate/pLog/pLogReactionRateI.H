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

#include "pLogReactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pLogReactionRate::pLogReactionRate
(
    const List<pCoeffs>& coeffs
)
:
    coeffs_(coeffs)
{}


inline Foam::pLogReactionRate::pLogReactionRate
(
    const speciesTable&,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    coeffs_(dict.lookup("coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::pLogReactionRate::preEvaluate() const
{}


inline void Foam::pLogReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::pLogReactionRate::ki
(
    const pCoeffs& coeffs,
    const scalar T
) const
{
    scalar k = coeffs.A;

    if (mag(coeffs.beta) > vSmall)
    {
        k *= pow(T, coeffs.beta);
    }

    if (mag(coeffs.Ta) > vSmall)
    {
        k *= exp(-coeffs.Ta/T);
    }

    return k;
}


inline Foam::scalar Foam::pLogReactionRate::lnki
(
    const pCoeffs& coeffs,
    const scalar T
) const
{
    scalar lnk = coeffs.lnA;

    if (mag(coeffs.beta) > vSmall)
    {
        lnk += coeffs.beta*log(T);
    }

    if (mag(coeffs.Ta) > vSmall)
    {
        lnk -= coeffs.Ta/T;
    }

    return lnk;
}


inline Foam::scalar Foam::pLogReactionRate::ddTi
(
    const pCoeffs& coeffs,
    const scalar T
) const
{
    return ki(coeffs, T)*(coeffs.beta + coeffs.Ta/T)/T;
}


inline Foam::scalar Foam::pLogReactionRate::fp
(
    const scalar p,
    const label pi
) const
{
    return
    (
        (log(p) - coeffs_[pi - 1].lnp)
       /(coeffs_[pi].lnp - coeffs_[pi - 1].lnp)
    );
}


inline Foam::scalar Foam::pLogReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    if (p <= coeffs_.first().p)
    {
        return ki(coeffs_.first(), T);
    }
    else if (p >= coeffs_.last().p)
    {
        return ki(coeffs_.last(), T);
    }
    else
    {
        label pi = 1;

        // Find the index of the pressure > p
        while (p > coeffs_[pi].p)
        {
            pi++;
        }

        // Calculate the log reaction rate coefficients
        // for the bounding pressures
        const scalar lnkim1(this->lnki(coeffs_[pi - 1], T));
        const scalar lnki(this->lnki(coeffs_[pi], T));

        // Interpolate the log reaction rate coefficients
        return exp(lnkim1 + fp(p, pi)*(lnki - lnkim1));
    }
}


inline Foam::scalar Foam::pLogReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    if (p <= coeffs_.first().p)
    {
        return ddTi(coeffs_.first(), T);
    }
    else if (p >= coeffs_.last().p)
    {
        return ddTi(coeffs_.last(), T);
    }
    else
    {
        label pi = 1;

        // Find the index of the pressure > p
        while (p > coeffs_[pi].p)
        {
            pi++;
        }

        // Calculate the log reaction rate coefficients
        // for the bounding pressures
        const scalar lnkim1(this->lnki(coeffs_[pi - 1], T));
        const scalar lnki(this->lnki(coeffs_[pi], T));

        // Calculate the interpolate factor
        const scalar fp = this->fp(p, pi);

        // Interpolate the log reaction rate coefficients
        const scalar k = exp(lnkim1 + fp*(lnki - lnkim1));

        // Multiply the interpolated log reaction rate coefficient
        // by the interpolated temperature derivative factor
        return k*
        (
            coeffs_[pi - 1].beta + fp*(coeffs_[pi].beta - coeffs_[pi - 1].beta)
          + (coeffs_[pi - 1].Ta + fp*(coeffs_[pi].Ta - coeffs_[pi - 1].Ta))/T
        )/T;
    }
}


inline bool Foam::pLogReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::pLogReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    NotImplemented;
}


inline void Foam::pLogReactionRate::write(Ostream& os) const
{
    writeEntry(os, "coeffs", coeffs_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pLogReactionRate& psrr
)
{
    psrr.write(os);
    return os;
}


inline bool Foam::operator==
(
    const pLogReactionRate::pCoeffs& ce1,
    const pLogReactionRate::pCoeffs& ce2
)
{
    return
    (
        ce1.p == ce2.p
     && ce1.A == ce2.A
     && ce1.beta == ce2.beta
     && ce1.Ta == ce2.Ta
    );
}


inline bool Foam::operator!=
(
    const pLogReactionRate::pCoeffs& ce1,
    const pLogReactionRate::pCoeffs& ce2
)
{
    return !(ce1 == ce2);
}


inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    pLogReactionRate::pCoeffs& ce
)
{
    is.readBegin("pCoeffs");
    is >> ce.p >> ce.A >> ce.beta >> ce.Ta;
    is.readEnd("pCoeffs");

    // Check state of Istream
    is.check("operator>>(Istream&, pLogReactionRate::pCoeffs&)");

    ce.lnp = log(ce.p);
    ce.lnA = log(ce.A);

    return is;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pLogReactionRate::pCoeffs& ce
)
{
    os  << token::BEGIN_LIST
        << ce.p << token::SPACE
        << ce.A << token::SPACE
        << ce.beta << token::SPACE
        << ce.Ta
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
