/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "inhomogeneousEGRMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::inhomogeneousEGRMixture<ThermoType>::inhomogeneousEGRMixture
(
    const dictionary& dict
)
:
    stoicRatio_("stoichiometricAirFuelMassRatio", dimless, dict),
    fuel_("fuel", dict.subDict("fuel")),
    oxidant_("oxidant", dict.subDict("oxidant")),
    products_("burntProducts", dict.subDict("burntProducts")),
    mixture_("mixture", fuel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::inhomogeneousEGRMixture<ThermoType>::fres
(
    const scalarFieldListSlice& Y
) const
{
    return max(Y[FT] - (1 - Y[FT] - Y[EGR])/stoicRatio_.value(), 0);
}


template<class ThermoType>
const ThermoType& Foam::inhomogeneousEGRMixture<ThermoType>::mixture
(
    const scalar ft,
    const scalar fu,
    const scalar egr
) const
{
    if (ft < 0.0001 && egr < 0.0001)
    {
        return oxidant_;
    }
    else
    {
        const scalar ox = 1 - ft - egr - (ft - fu)*stoicRatio_.value();
        const scalar pr = 1 - fu - ox;

        mixture_ = fu*fuel_;
        mixture_ += ox*oxidant_;
        mixture_ += pr*products_;

        return mixture_;
    }
}


template<class ThermoType>
const typename Foam::inhomogeneousEGRMixture<ThermoType>::thermoMixtureType&
Foam::inhomogeneousEGRMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], Y[FU], Y[EGR]);
}


template<class ThermoType>
const typename Foam::inhomogeneousEGRMixture<ThermoType>::transportMixtureType&
Foam::inhomogeneousEGRMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], Y[FU], Y[EGR]);
}


template<class ThermoType>
const typename Foam::inhomogeneousEGRMixture<ThermoType>::transportMixtureType&
Foam::inhomogeneousEGRMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


template<class ThermoType>
const typename Foam::inhomogeneousEGRMixture<ThermoType>::thermoType&
Foam::inhomogeneousEGRMixture<ThermoType>::reactants
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], Y[FT], Y[EGR]);
}


template<class ThermoType>
const typename Foam::inhomogeneousEGRMixture<ThermoType>::thermoType&
Foam::inhomogeneousEGRMixture<ThermoType>::products
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], fres(Y), Y[EGR]);
}


template<class ThermoType>
void Foam::inhomogeneousEGRMixture<ThermoType>::reset
(
    PtrList<volScalarField>& Y
) const
{
    const volScalarField& fu = Y[FU];
    volScalarField& ft = Y[FT];
    volScalarField& b = Y[B];
    volScalarField& egr = Y[EGR];

    for (label i=0; i<=fu.nOldTimes(); i++)
    {
        egr.oldTimeRef(i) += (ft.oldTime(i) - fu.oldTime(i))*(1 + stoicRatio_);
        ft.oldTimeRef(i) = fu.oldTime(i);
        b.oldTimeRef(i) = 1;
    }
}


template<class ThermoType>
void Foam::inhomogeneousEGRMixture<ThermoType>::read(const dictionary& dict)
{
    stoicRatio_ =
        dimensionedScalar("stoichiometricAirFuelMassRatio", dimless, dict);

    fuel_ = ThermoType("fuel", dict.subDict("fuel"));
    oxidant_ = ThermoType("oxidant", dict.subDict("oxidant"));
    products_ = ThermoType("burntProducts", dict.subDict("burntProducts"));
}


// ************************************************************************* //
