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

#include "ThirdBody_reactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Rate>
inline Foam::reactionRates::ThirdBody<Rate>::ThirdBody
(
    const Rate& k,
    const thirdBodyEfficiencies& tbes
)
:
    k_(k),
    thirdBodyEfficiencies_(tbes)
{}


template<class Rate>
inline Foam::reactionRates::ThirdBody<Rate>::ThirdBody
(
    const speciesTable& species,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k_(species, dims/(dimensions::moles/dimensions::volume), dict),
    thirdBodyEfficiencies_(species, dict)
{}


template<class Rate>
inline Foam::reactionRates::ThirdBody<Rate>::ThirdBody
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    k_(species, ob, dims/(dimensions::moles/dimensions::volume), dict),
    thirdBodyEfficiencies_(species, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Rate>
inline void Foam::reactionRates::ThirdBody<Rate>::preEvaluate() const
{
    k_.preEvaluate();
}


template<class Rate>
inline void Foam::reactionRates::ThirdBody<Rate>::postEvaluate() const
{
    k_.postEvaluate();
}


template<class Rate>
inline Foam::scalar Foam::reactionRates::ThirdBody<Rate>::k
(
    const reaction::state& s
) const
{
    return thirdBodyEfficiencies_.M(s)*k_.k(s);
}


template<class Rate>
inline Foam::scalar Foam::reactionRates::ThirdBody<Rate>::ddT
(
    const reaction::state& s
) const
{
    return thirdBodyEfficiencies_.M(s)*k_.ddT(s);
}


template<class Rate>
inline void Foam::reactionRates::ThirdBody<Rate>::ddc
(
    const reaction::state& s,
    scalarField& ddc
) const
{
    ddc = thirdBodyEfficiencies_.dMdc();
    ddc *= k_.k(s);
}


template<class Rate>
inline void Foam::reactionRates::ThirdBody<Rate>::write(Ostream& os) const
{
    k_.write(os);
    thirdBodyEfficiencies_.write(os);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Rate>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactionRates::ThirdBody<Rate>& rr
)
{
    rr.write(os);
    return os;
}


// ************************************************************************* //
