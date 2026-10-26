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

#include "ChemicallyActivated_reactionRate.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate, class ChemicalActivationFunction>
inline Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::ChemicallyActivated
(
    const Rate& k0,
    const Rate& kInf,
    const ChemicalActivationFunction& F,
    const thirdBodyEfficiencies& tbes
)
:
    k0_(k0),
    kInf_(kInf),
    F_(F),
    thirdBodyEfficiencies_(tbes)
{}


template<class Rate, class ChemicalActivationFunction>
inline Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::ChemicallyActivated
(
    const speciesTable& species,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k0_(species, dims, dict.subDict("k0")),
    kInf_
    (
        species,
        dims*(dimensions::moles/dimensions::volume),
        dict.subDict("kInf")
    ),
    F_(dict.subDict("F")),
    thirdBodyEfficiencies_(species, dict.subDict("thirdBodyEfficiencies"))
{}


template<class Rate, class ChemicalActivationFunction>
inline Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::ChemicallyActivated
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k0_(species, ob, dims, dict.subDict("k0")),
    kInf_
    (
        species,
        ob,
        dims*(dimensions::moles/dimensions::volume),
        dict.subDict("kInf")
    ),
    F_(dict.subDict("F")),
    thirdBodyEfficiencies_(species, dict.subDict("thirdBodyEfficiencies"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate, class ChemicalActivationFunction>
inline void Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::preEvaluate() const
{
    k0_.preEvaluate();
    kInf_.preEvaluate();
}


template<class Rate, class ChemicalActivationFunction>
inline void Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::postEvaluate() const
{
    k0_.postEvaluate();
    kInf_.postEvaluate();
}


template<class Rate, class ChemicalActivationFunction>
inline Foam::scalar Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::k(const reaction::state& s) const
{
    const scalar k0 = k0_.k(s);
    const scalar kInf = kInf_.k(s);
    const scalar M = thirdBodyEfficiencies_.M(s);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_.F(s.T, Pr);

    return k0/(1 + Pr)*F;
}


template<class Rate, class ChemicalActivationFunction>
inline Foam::scalar Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::ddT(const reaction::state& s) const
{
    const scalar k0 = k0_.k(s);
    const scalar kInf = kInf_.k(s);
    const scalar M = thirdBodyEfficiencies_.M(s);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_.F(s.T, Pr);

    const scalar dk0dT = k0_.ddT(s);

    #if CHEMICAL_ACTIVATION_FUNCTION_JACOBIAN

    const scalar dkInfdT = kInf_.ddT(s);
    const scalar dPrdT = (M*dk0dT - Pr*dkInfdT)/kInf;
    const scalar dFdT = F_.ddT(s.T, Pr, F) + F_.ddPr(s.T, Pr, F)*dPrdT;

    return
        dk0dT/(1 + Pr)*F
      - k0*dPrdT/sqr(1 + Pr)*F
      + k0/(1 + Pr)*dFdT;

    #else

    return dk0dT/(1 + Pr)*F;

    #endif
}


template<class Rate, class ChemicalActivationFunction>
inline void Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::ddc
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

    #if CHEMICAL_ACTIVATION_FUNCTION_JACOBIAN

    scalarField dk0dc(c.size(), 0);
    k0_.ddc(s, dk0dc);
    scalarField dkInfdc(c.size(), 0);
    kInf_.ddc(s, dkInfdc);
    tmp<scalarField> tdMdc = thirdBodyEfficiencies_.dMdc(s);
    const scalarField& dMdc = tdMdc();
    scalarField dPrdc((M*dk0dc - Pr*dkInfdc + k0*dMdc)/kInf);
    const scalar dFdPr = F_.ddPr(s.T, Pr, F);

    ddc =
        dk0dc/(1 + Pr)*F
      - k0*dPrdc/sqr(1 + Pr)*F
      + k0/(1 + Pr)*dFdPr*dPrdc;

    #else

    k0_.ddc(s, ddc);
    ddc *= 1/(1 + Pr)*F;

    #endif
}


template<class Rate, class ChemicalActivationFunction>
inline void Foam::reactionRates::ChemicallyActivated
<
    Rate,
    ChemicalActivationFunction
>::write(Ostream& os) const
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

template<class Rate, class ChemicalActivationFunction>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::ChemicallyActivated
        <Rate, ChemicalActivationFunction>& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
