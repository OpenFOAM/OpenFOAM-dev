/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "FallOff_reactionRate.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define FALL_OFF_FUNCTION_JACOBIAN false

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate, class FallOffFunction>
inline Foam::reactionRates::FallOff<Rate, FallOffFunction>::FallOff
(
    const Rate& k0,
    const Rate& kInf,
    const FallOffFunction& F,
    const thirdBodyEfficiencies& tbes
)
:
    k0_(k0),
    kInf_(kInf),
    F_(F),
    thirdBodyEfficiencies_(tbes)
{}


template<class Rate, class FallOffFunction>
inline Foam::reactionRates::FallOff<Rate, FallOffFunction>::FallOff
(
    const speciesTable& species,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k0_
    (
        species,
        dims/(dimensions::moles/dimensions::volume),
        dict.subDict("k0")
    ),
    kInf_(species, dims, dict.subDict("kInf")),
    F_(dict.subDict("F")),
    thirdBodyEfficiencies_(species, dict.subDict("thirdBodyEfficiencies"))
{}


template<class Rate, class FallOffFunction>
inline Foam::reactionRates::FallOff<Rate, FallOffFunction>::FallOff
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k0_
    (
        species,
        ob,
        dims/(dimensions::moles/dimensions::volume),
        dict.subDict("k0")
    ),
    kInf_(species, ob, dims, dict.subDict("kInf")),
    F_(dict.subDict("F")),
    thirdBodyEfficiencies_(species, dict.subDict("thirdBodyEfficiencies"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate, class FallOffFunction>
inline void
Foam::reactionRates::FallOff<Rate, FallOffFunction>::preEvaluate() const
{
    k0_.preEvaluate();
    kInf_.preEvaluate();
}


template<class Rate, class FallOffFunction>
inline void
Foam::reactionRates::FallOff<Rate, FallOffFunction>::postEvaluate() const
{
    k0_.postEvaluate();
    kInf_.postEvaluate();
}


template<class Rate, class FallOffFunction>
inline Foam::scalar Foam::reactionRates::FallOff<Rate, FallOffFunction>::k
(
    const reaction::state& s
) const
{
    const scalar k0 = k0_.k(s);
    const scalar kInf = kInf_.k(s);
    const scalar M = thirdBodyEfficiencies_.M(s);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_.F(s.T, Pr);

    return kInf*(Pr/(1 + Pr))*F;
}


template<class Rate, class FallOffFunction>
inline Foam::scalar Foam::reactionRates::FallOff<Rate, FallOffFunction>::ddT
(
    const reaction::state& s
) const
{
    const scalar k0 = k0_.k(s);
    const scalar kInf = kInf_.k(s);
    const scalar M = thirdBodyEfficiencies_.M(s);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_.F(s.T, Pr);

    const scalar dkInfdT = kInf_.ddT(s);

    #if FALL_OFF_FUNCTION_JACOBIAN

    const scalar dk0dT = k0_.ddT(s);
    const scalar dPrdT = (M*dk0dT - Pr*dkInfdT)/kInf;
    const scalar dFdT = F_.ddT(s.T, Pr, F) + F_.ddPr(s.T, Pr, F)*dPrdT;

    return
        dkInfdT*(Pr/(1 + Pr))*F
      + kInf*dPrdT/sqr(1 + Pr)*F
      + kInf*(Pr/(1 + Pr))*dFdT;

    #else

    return dkInfdT*(Pr/(1 + Pr))*F;

    #endif
}


template<class Rate, class FallOffFunction>
inline void Foam::reactionRates::FallOff<Rate, FallOffFunction>::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    const scalar k0 = k0_.k(s);
    const scalar kInf = kInf_.k(s);
    const scalar M = thirdBodyEfficiencies_.M(s);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_.F(s.T, Pr);

    #if FALL_OFF_FUNCTION_JACOBIAN

    scalarField dk0dc(s.c.size(), 0);
    k0_.ddc(s, dk0dc);
    scalarField dkInfdc(s.c.size(), 0);
    kInf_.ddc(s, dkInfdc);
    tmp<scalarField> tdMdc = thirdBodyEfficiencies_.dMdc();
    const scalarField& dMdc = tdMdc();
    scalarField dPrdc((M*dk0dc - Pr*dkInfdc + k0*dMdc)/kInf);
    const scalar dFdPr = F_.ddPr(s.T, Pr, F);

    ddc =
        dkInfdc*(Pr/(1 + Pr))*F
      + kInf*dPrdc/sqr(1 + Pr)*F
      + kInf*(Pr/(1 + Pr))*dFdPr*dPrdc;

    #else

    kInf_.ddc(s, ddc);
    ddc *= (Pr/(1 + Pr))*F;

    #endif
}


template<class Rate, class FallOffFunction>
inline void Foam::reactionRates::FallOff<Rate, FallOffFunction>::write
(
    Ostream& os
) const
{
    {
        const delimitDictionary delimitK0(os, "k0");
        k0_.write(os);
    }
    {
        const delimitDictionary delimitKInf(os, "kInf");
        kInf_.write(os);
    }
    {
        const delimitDictionary delimitF(os, "F");
        F_.write(os);
    }
    {
        const delimitDictionary delimitTbe(os, "thirdBodyEfficiencies");
        thirdBodyEfficiencies_.write(os);
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Rate, class FallOffFunction>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::FallOff<Rate, FallOffFunction>& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
