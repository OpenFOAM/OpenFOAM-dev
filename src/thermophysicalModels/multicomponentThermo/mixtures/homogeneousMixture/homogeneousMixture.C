/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "homogeneousMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::homogeneousMixture<ThermoType>::homogeneousMixture
(
    const dictionary& dict
)
:
    reactants_("reactants", dict.subDict("reactants")),
    products_("products", dict.subDict("products")),
    mixture_("mixture", reactants_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::homogeneousMixture<ThermoType>::mixture
(
    const scalar b
) const
{
    if (b > 0.999)
    {
        return reactants_;
    }
    else if (b < 0.001)
    {
        return products_;
    }
    else
    {
        mixture_ = b*reactants_;
        mixture_ += (1 - b)*products_;

        return mixture_;
    }
}


template<class ThermoType>
const typename Foam::homogeneousMixture<ThermoType>::thermoMixtureType&
Foam::homogeneousMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[B]);
}


template<class ThermoType>
const typename Foam::homogeneousMixture<ThermoType>::transportMixtureType&
Foam::homogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[B]);
}


template<class ThermoType>
const typename Foam::homogeneousMixture<ThermoType>::transportMixtureType&
Foam::homogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


template<class ThermoType>
const typename Foam::homogeneousMixture<ThermoType>::thermoType&
Foam::homogeneousMixture<ThermoType>::reactants
(
    const scalarFieldListSlice& Y
) const
{
    return reactants_;
}


template<class ThermoType>
const typename Foam::homogeneousMixture<ThermoType>::thermoType&
Foam::homogeneousMixture<ThermoType>::products
(
    const scalarFieldListSlice& Y
) const
{
    return products_;
}


template<class ThermoType>
void Foam::homogeneousMixture<ThermoType>::read(const dictionary& dict)
{
    reactants_ = ThermoType("reactants", dict.subDict("reactants"));
    products_ = ThermoType("products", dict.subDict("products"));
}


// ************************************************************************* //
