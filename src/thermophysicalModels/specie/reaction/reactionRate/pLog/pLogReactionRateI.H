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
    const reaction::state& s
) const
{
    return coeffs.A*exp(coeffs.beta*s.lnT - coeffs.Ta*s.rT);
}


inline Foam::scalar Foam::pLogReactionRate::lnki
(
    const pCoeffs& coeffs,
    const reaction::state& s
) const
{
    return coeffs.lnA + coeffs.beta*s.lnT - coeffs.Ta*s.rT;
}


inline Foam::scalar Foam::pLogReactionRate::ddTi
(
    const pCoeffs& coeffs,
    const reaction::state& s
) const
{
    return ki(coeffs, s)*(coeffs.beta + coeffs.Ta*s.rT)*s.rT;
}


inline Foam::scalar Foam::pLogReactionRate::fp
(
    const scalar lnp,
    const label pi
) const
{
    return
    (
        (lnp - coeffs_[pi - 1].lnp)
       /(coeffs_[pi].lnp - coeffs_[pi - 1].lnp)
    );
}


inline Foam::scalar Foam::pLogReactionRate::k
(
    const reaction::state& s
) const
{
    if (s.p <= coeffs_.first().p)
    {
        return ki(coeffs_.first(), s);
    }
    else if (s.p >= coeffs_.last().p)
    {
        return ki(coeffs_.last(), s);
    }
    else
    {
        label pi = 1;

        // Find the index of the pressure > p
        while (s.p > coeffs_[pi].p)
        {
            pi++;
        }

        // Calculate the log reaction rate coefficients
        // for the bounding pressures
        const scalar lnkim1(this->lnki(coeffs_[pi - 1], s));
        const scalar lnki(this->lnki(coeffs_[pi], s));

        // Interpolate the log reaction rate coefficients
        return exp(lnkim1 + fp(s.lnp, pi)*(lnki - lnkim1));
    }
}


inline Foam::scalar Foam::pLogReactionRate::ddT
(
    const reaction::state& s
) const
{
    if (s.p <= coeffs_.first().p)
    {
        return ddTi(coeffs_.first(), s);
    }
    else if (s.p >= coeffs_.last().p)
    {
        return ddTi(coeffs_.last(), s);
    }
    else
    {
        label pi = 1;

        // Find the index of the pressure > p
        while (s.p > coeffs_[pi].p)
        {
            pi++;
        }

        // Calculate the log reaction rate coefficients
        // for the bounding pressures
        const scalar lnkim1(this->lnki(coeffs_[pi - 1], s));
        const scalar lnki(this->lnki(coeffs_[pi], s));

        // Calculate the interpolate factor
        const scalar fp = this->fp(s.lnp, pi);

        // Interpolate the log reaction rate coefficients
        const scalar k = exp(lnkim1 + fp*(lnki - lnkim1));

        // Multiply the interpolated log reaction rate coefficient
        // by the interpolated temperature derivative factor
        return k*
        (
            coeffs_[pi - 1].beta + fp*(coeffs_[pi].beta - coeffs_[pi - 1].beta)
          + (coeffs_[pi - 1].Ta + fp*(coeffs_[pi].Ta - coeffs_[pi - 1].Ta))*s.rT
        )*s.rT;
    }
}


inline bool Foam::pLogReactionRate::hasDdc() const
{
    return false;
}


inline void Foam::pLogReactionRate::ddc
(
    const reaction::state& s,
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
