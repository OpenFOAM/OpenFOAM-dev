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

#include "uInhomogeneousEGRMixture.H"
#include "bInhomogeneousMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::uInhomogeneousEGRMixture<ThermoType>::uInhomogeneousEGRMixture
(
    const dictionary& dict
)
:
    species_({"fu", "egr"}),
    stoicRatio_(dict.lookup<scalar>("stoichiometricAirFuelMassRatio")),
    fuel_("fuel", dict.subDict("fuel")),
    oxidant_("oxidant", dict.subDict("oxidant")),
    products_("products", dict.subDict("products")),
    active_(2, true),
    mixture_("mixture", fuel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::uInhomogeneousEGRMixture<ThermoType>::specieThermo
(
    const label speciei
) const
{
    switch (speciei)
    {
        case FU:
            return fuel_;
            break;

        case EGR:
            return products_;
            break;

        default:
            FatalErrorInFunction
                << "Cannot return specieThermo for specie " << speciei
                << exit(FatalError);
            return products_;
            break;
    }
}


template<class ThermoType>
Foam::scalar Foam::uInhomogeneousEGRMixture<ThermoType>::Phi
(
    const scalarFieldListSlice& Yu
) const
{
    const scalar ft = Yu[FU] + Yu[EGR]/(stoicRatio_ + 1);
    return stoicRatio_*ft/max(1 - ft, small);
}


template<class ThermoType>
const ThermoType& Foam::uInhomogeneousEGRMixture<ThermoType>::mixture
(
    const scalar fu,
    const scalar egr
) const
{
    if (fu < 0.0001 && egr < 0.0001)
    {
        return oxidant_;
    }
    else
    {
        const scalar ox = 1 - fu - egr;

        mixture_ = fu*fuel_;
        mixture_ += ox*oxidant_;
        mixture_ += egr*products_;

        return mixture_;
    }
}


template<class ThermoType>
const typename Foam::uInhomogeneousEGRMixture<ThermoType>::thermoMixtureType&
Foam::uInhomogeneousEGRMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FU], Y[EGR]);
}


template<class ThermoType>
const typename Foam::uInhomogeneousEGRMixture<ThermoType>::transportMixtureType&
Foam::uInhomogeneousEGRMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FU], Y[EGR]);
}


template<class ThermoType>
const typename Foam::uInhomogeneousEGRMixture<ThermoType>::transportMixtureType&
Foam::uInhomogeneousEGRMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


template<class ThermoType>
Foam::PtrList<Foam::volScalarField::Internal>
Foam::uInhomogeneousEGRMixture<ThermoType>::prompt
(
    const PtrList<volScalarField>& Yu
) const
{
    PtrList<volScalarField::Internal> Yp(1);
    Yp.set(bInhomogeneousMixture<ThermoType>::FT, Yu[FU]());

    return Yp;
}


template<class ThermoType>
void Foam::uInhomogeneousEGRMixture<ThermoType>::reset
(
    const volScalarField& b,
    PtrList<volScalarField>& Yu,
    const volScalarField& c,
    const PtrList<volScalarField>& Yb
) const
{
    volScalarField& fuu = Yu[FU];
    volScalarField& egr = Yu[EGR];

    const volScalarField& ftb = Yb[bInhomogeneousMixture<ThermoType>::FT];

    for (label i=0; i<=fuu.nOldTimes(); i++)
    {
        const volScalarField fub
        (
            max
            (
                ftb.oldTime(i) - (scalar(1) - ftb.oldTime(i))/stoicRatio_,
                scalar(0)
            )
        );
        const volScalarField oxb
        (
            1 - ftb.oldTime(i) - (ftb.oldTime(i) - fub)*stoicRatio_
        );

        fuu.oldTimeRef(i) = b.oldTime(i)*fuu.oldTime(i) + c.oldTime(i)*fub;
        egr.oldTimeRef(i) =
            b.oldTime(i)*egr.oldTime(i) + c.oldTime(i)*(1 - fub - oxb);
    }
}


template<class ThermoType>
void Foam::uInhomogeneousEGRMixture<ThermoType>::read(const dictionary& dict)
{
    stoicRatio_ = dict.lookup<scalar>("stoichiometricAirFuelMassRatio");
    fuel_ = ThermoType("fuel", dict.subDict("fuel"));
    oxidant_ = ThermoType("oxidant", dict.subDict("oxidant"));
    products_ = ThermoType("products", dict.subDict("products"));
}


// ************************************************************************* //
