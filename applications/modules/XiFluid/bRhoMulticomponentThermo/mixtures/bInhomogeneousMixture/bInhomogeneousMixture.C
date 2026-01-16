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

#include "bInhomogeneousMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::bInhomogeneousMixture<ThermoType>::bInhomogeneousMixture
(
    const dictionary& dict
)
:
    stoicRatio_(dict.lookup<scalar>("stoichiometricAirFuelMassRatio")),
    fuel_("fuel", dict.subDict("fuel")),
    oxidant_("oxidant", dict.subDict("oxidant")),
    products_("products", dict.subDict("products")),
    active_(1, true),
    mixture_("mixture", fuel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::bInhomogeneousMixture<ThermoType>::mixture
(
    const scalar ft
) const
{
    if (ft < 0.0001)
    {
        return oxidant_;
    }
    else
    {
        const scalar fu = max(ft - (scalar(1) - ft)/stoicRatio_, scalar(0));
        const scalar ox = 1 - ft - (ft - fu)*stoicRatio_;
        const scalar pr = 1 - fu - ox;

        mixture_ = fu*fuel_;
        mixture_ += ox*oxidant_;
        mixture_ += pr*products_;

        return mixture_;
    }
}


template<class ThermoType>
const typename Foam::bInhomogeneousMixture<ThermoType>::thermoMixtureType&
Foam::bInhomogeneousMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT]);
}


template<class ThermoType>
const typename Foam::bInhomogeneousMixture<ThermoType>::transportMixtureType&
Foam::bInhomogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT]);
}


template<class ThermoType>
const typename Foam::bInhomogeneousMixture<ThermoType>::transportMixtureType&
Foam::bInhomogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


template<class ThermoType>
void Foam::bInhomogeneousMixture<ThermoType>::read(const dictionary& dict)
{
    stoicRatio_ = dict.lookup<scalar>("stoichiometricAirFuelMassRatio");
    fuel_ = ThermoType("fuel", dict.subDict("fuel"));
    oxidant_ = ThermoType("oxidant", dict.subDict("oxidant"));
    products_ = ThermoType("products", dict.subDict("products"));
}


// ************************************************************************* //
