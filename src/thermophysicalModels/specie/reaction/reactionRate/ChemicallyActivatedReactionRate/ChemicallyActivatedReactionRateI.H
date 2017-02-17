/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    k0_(species, dict),
    kInf_(species, dict),
    F_(dict),
    thirdBodyEfficiencies_(species, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionRate, class ChemicallyActivationFunction>
inline Foam::scalar Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    scalar k0 = k0_(p, T, c);
    scalar kInf = kInf_(p, T, c);

    scalar Pr = k0*thirdBodyEfficiencies_.M(c)/kInf;

    return k0*(1/(1 + Pr))*F_(T, Pr);
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline void Foam::ChemicallyActivatedReactionRate
<
    ReactionRate,
    ChemicallyActivationFunction
>::write(Ostream& os) const
{
    k0_.write(os);
    kInf_.write(os);
    F_.write(os);
    thirdBodyEfficiencies_.write(os);
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
