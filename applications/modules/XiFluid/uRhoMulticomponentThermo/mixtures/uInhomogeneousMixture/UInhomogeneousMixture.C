/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "UInhomogeneousMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::UInhomogeneousMixture<ThermoType>::UInhomogeneousMixture
(
    const dictionary& dict
)
:
    uInhomogeneousMixture(dict),
    fuel_("fuel", dict.subDict("fuel")),
    oxidant_("oxidant", dict.subDict("oxidant")),
    mixture_("mixture", fuel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::UInhomogeneousMixture<ThermoType>::mixture
(
    const scalar fu
) const
{
    if (fu < 0.0001)
    {
        return oxidant_;
    }
    else
    {
        const scalar ox = 1 - fu;

        mixture_ = fu*fuel_;
        mixture_ += ox*oxidant_;

        return mixture_;
    }
}


template<class ThermoType>
const typename Foam::UInhomogeneousMixture<ThermoType>::thermoMixtureType&
Foam::UInhomogeneousMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FU]);
}


template<class ThermoType>
const typename Foam::UInhomogeneousMixture<ThermoType>::transportMixtureType&
Foam::UInhomogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FU]);
}


template<class ThermoType>
const typename Foam::UInhomogeneousMixture<ThermoType>::transportMixtureType&
Foam::UInhomogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


template<class ThermoType>
void Foam::UInhomogeneousMixture<ThermoType>::read(const dictionary& dict)
{
    uInhomogeneousMixture::read(dict);
    fuel_ = ThermoType("fuel", dict.subDict("fuel"));
    oxidant_ = ThermoType("oxidant", dict.subDict("oxidant"));
}


// ************************************************************************* //
