/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "reaction.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline Foam::dimensionSet Foam::reaction::kDimensions
(
    const List<specieCoeffs>& coeffs
)
{
    scalar order = 0;
    forAll(coeffs, i)
    {
        order += coeffs[i].exponent;
    }

    return
        pow(dimensions::moles/dimensions::volume, 1 - order)
       /dimensions::time;
}


inline Foam::scalar Foam::reaction::cPowExponentProduct
(
    const List<specieCoeffs>& coeffs,
    const state& s
)
{
    scalar result = pow(s.c[coeffs[0].index], coeffs[0].exponent);

    for(label i = 1; i < coeffs.size(); ++ i)
    {
        result *= pow(s.c[coeffs[i].index], coeffs[i].exponent);
    }

    return result;
}


inline Foam::scalar Foam::reaction::cPowExponentProduct
(
    const List<specieCoeffs>& coeffs,
    const state& s,
    scalarField& dResultDc
)
{
    const label n = coeffs.size();

    const scalar cPowExponent =
        pow(s.c[coeffs[0].index], coeffs[0].exponent);

    scalar result = cPowExponent;
    dResultDc[0] = dpowdc(s.c[coeffs[0].index], coeffs[0].exponent);
    for (label j = 1; j < n; ++ j) dResultDc[j] = cPowExponent;

    for (label i = 1; i < n; ++ i)
    {
        const scalar cPowExponent =
            pow(s.c[coeffs[i].index], coeffs[i].exponent);

        result *= cPowExponent;
        for (label j = 0; j < i; ++ j) dResultDc[j] *= cPowExponent;
        dResultDc[i] *= dpowdc(s.c[coeffs[i].index], coeffs[i].exponent);
        for (label j = i + 1; j < n; ++ j) dResultDc[j] *= cPowExponent;
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::reaction::omega(const state& s) const
{
    scalar omegaf, omegar;

    return omega(s, omegaf, omegar);
}


template<class C2S>
inline void Foam::reaction::dNdtByV
(
    const scalar omega,
    scalarField& dNdtByV,
    const C2S& c2s
) const
{
    forAll(lhs(), i)
    {
        dNdtByV[c2s[lhs()[i].index]] -= lhs()[i].stoichCoeff*omega;
    }
    forAll(rhs(), i)
    {
        dNdtByV[c2s[rhs()[i].index]] += rhs()[i].stoichCoeff*omega;
    }
}


inline void Foam::reaction::dNdtByV
(
    const scalar omega,
    scalarField& dNdtByV,
    const bool reduced,
    const List<label>& c2s
) const
{
    if (reduced)
    {
        this->dNdtByV(omega, dNdtByV, c2s);
    }
    else
    {
        this->dNdtByV(omega, dNdtByV, noc2s());
    }
}


template<class C2S>
inline void Foam::reaction::dNdtByV
(
    const state& s,
    scalarField& dNdtByV,
    const C2S& c2s
) const
{
    return this->dNdtByV(reaction::omega(s), dNdtByV, c2s);
}


inline void Foam::reaction::dNdtByV
(
    const state& s,
    scalarField& dNdtByV,
    const bool reduced,
    const List<label>& c2s
) const
{
    if (reduced)
    {
        this->dNdtByV(s, dNdtByV, c2s);
    }
    else
    {
        this->dNdtByV(s, dNdtByV, noc2s());
    }
}


template<class C2S>
void Foam::reaction::ddNdtByVdcT
(
    const state& s,
    scalarField& dNdtByV,
    scalarSquareMatrix& ddNdtByVdcT,
    const C2S& c2s,
    const label Tj
) const
{
    // Net rates and derivatives
    scalar dOmegaDT;
    scalarField& dOmegaDc = s.specieWork_.elmt<0>();
    const scalar omega = this->omega(s, dOmegaDT, dOmegaDc);

    // Rates
    this->dNdtByV(omega, dNdtByV, c2s);

    // Derivatives w.r.t temperature
    forAll(lhs(), i)
    {
        ddNdtByVdcT(c2s[lhs()[i].index], Tj) -=
            lhs()[i].stoichCoeff*dOmegaDT;
    }
    forAll(rhs(), i)
    {
        ddNdtByVdcT(c2s[rhs()[i].index], Tj) +=
            rhs()[i].stoichCoeff*dOmegaDT;
    }

    // Derivatives w.r.t concentration
    auto addDdc = [&](const label cj)
    {
        forAll(lhs(), i)
        {
            ddNdtByVdcT(c2s[lhs()[i].index], c2s[cj]) -=
                lhs()[i].stoichCoeff*dOmegaDc[cj];
        }
        forAll(rhs(), i)
        {
            ddNdtByVdcT(c2s[rhs()[i].index], c2s[cj]) +=
                rhs()[i].stoichCoeff*dOmegaDc[cj];
        }
    };
    switch (density())
    {
        case jacobianDensity::reactants:
            forAll(lhs(), j) addDdc(lhs()[j].index);
            break;
        case jacobianDensity::reactantsAndProducts:
            forAll(specieIndices(), j) addDdc(specieIndices()[j]);
            break;
        case jacobianDensity::all:
            forAll(s.c, j) if (c2s[j] != -1) addDdc(j);
            break;
    }
}


void Foam::reaction::ddNdtByVdcT
(
    const state& s,
    scalarField& dNdtByV,
    scalarSquareMatrix& ddNdtByVdcT,
    const bool reduced,
    const List<label>& c2s,
    const label Tj
) const
{
    if (reduced)
    {
        this->ddNdtByVdcT(s, dNdtByV, ddNdtByVdcT, c2s, Tj);
    }
    else
    {
        this->ddNdtByVdcT(s, dNdtByV, ddNdtByVdcT, noc2s(), Tj);
    }
}


inline void Foam::reaction::state::reset
(
    const scalar p_,
    const scalar T_,
    const label li_
)
{
    p = p_;
    T = T_;
    li = li_;
    lnp = log(p);
    rT = 1/T;
    lnT = log(T);
    C = sum(c);
}


inline void Foam::reaction::state::reset
(
    const scalar p_,
    const scalar T_,
    const scalarField& c_,
    const label li_
)
{
    c = c_;
    reset(p_, T_, li_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::reaction::noc2s::operator[](const label c) const
{
    return c;
}


// ************************************************************************* //
