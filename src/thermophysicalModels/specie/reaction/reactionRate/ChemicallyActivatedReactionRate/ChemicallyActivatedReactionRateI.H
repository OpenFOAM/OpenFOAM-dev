/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#define CHEMICALLY_ACTIVATION_FUNCTION_JACOBIAN false

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ReactionRate, class ChemicallyActivationFunction>
inline Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::ChemicallyActivatedReactionRate
(
    const ReactionRate& k0,
    const ReactionRate& kInf,
    const ChemicallyActivationFunction& F,
    const thirdBodyEfficiencies& tbes
)
:
    k0_(k0),
    kInf_(kInf),
    F_(F),
    thirdBodyEfficiencies_(tbes)
{}


template<class ReactionRate, class ChemicallyActivationFunction>
inline Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::ChemicallyActivatedReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    k0_(species, dict.subDict("k0")),
    kInf_(species, dict.subDict("kInf")),
    F_(dict.subDict("F")),
    thirdBodyEfficiencies_(species, dict.subDict("thirdBodyEfficiencies"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionRate, class ChemicallyActivationFunction>
inline void Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::preEvaluate() const
{
    k0_.preEvaluate();
    kInf_.preEvaluate();
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline void Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::postEvaluate() const
{
    k0_.postEvaluate();
    kInf_.postEvaluate();
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline Foam::scalar Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    const scalar k0 = k0_(p, T, c, li);
    const scalar kInf = kInf_(p, T, c, li);
    const scalar M = thirdBodyEfficiencies_.M(c);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_(T, Pr);

    return k0/(1 + Pr)*F;
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline Foam::scalar Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    const scalar k0 = k0_(p, T, c, li);
    const scalar kInf = kInf_(p, T, c, li);
    const scalar M = thirdBodyEfficiencies_.M(c);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_(T, Pr);

    const scalar dk0dT = k0_.ddT(p, T, c, li);

    #if CHEMICALLY_ACTIVATION_FUNCTION_JACOBIAN

    const scalar dkInfdT = kInf_.ddT(p, T, c, li);
    const scalar dPrdT = (M*dk0dT - Pr*dkInfdT)/kInf;
    const scalar dFdT = F_.ddT(T, Pr, F) + F_.ddPr(T, Pr, F)*dPrdT;

    return
        dk0dT/(1 + Pr)*F
      - k0*dPrdT/sqr(1 + Pr)*F
      + k0/(1 + Pr)*dFdT;

    #else

    return dk0dT/(1 + Pr)*F;

    #endif
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline bool Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::hasDdc() const
{
    return true;
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline void Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    const scalar k0 = k0_(p, T, c, li);
    const scalar kInf = kInf_(p, T, c, li);
    const scalar M = thirdBodyEfficiencies_.M(c);
    const scalar Pr = k0/kInf*M;
    const scalar F = F_(T, Pr);

    #if CHEMICALLY_ACTIVATION_FUNCTION_JACOBIAN

    scalarField dk0dc(c.size(), 0);
    k0_.ddc(p, T, c, li, dk0dc);
    scalarField dkInfdc(c.size(), 0);
    kInf_.ddc(p, T, c, li, dkInfdc);
    tmp<scalarField> tdMdc = thirdBodyEfficiencies_.dMdc(c);
    const scalarField& dMdc = tdMdc();
    scalarField dPrdc((M*dk0dc - Pr*dkInfdc + k0*dMdc)/kInf);
    const scalar dFdPr = F_.ddPr(T, Pr, F);

    ddc =
        dk0dc/(1 + Pr)*F
      - k0*dPrdc/sqr(1 + Pr)*F
      + k0/(1 + Pr)*dFdPr*dPrdc;

    #else

    k0_.ddc(p, T, c, li, ddc);

    ddc *= 1/(1 + Pr)*F;

    #endif
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline void Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::write(Ostream& os) const
{
    os  << indent << "k0" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    k0_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "kInf" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    kInf_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "F" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    F_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "thirdBodyEfficiencies" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    thirdBodyEfficiencies_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ChemicallyActivatedReactionRate
        <ReactionRate, ChemicallyActivationFunction>& carr
)
{
    carr.write(os);
    return os;
}


// ************************************************************************* //
