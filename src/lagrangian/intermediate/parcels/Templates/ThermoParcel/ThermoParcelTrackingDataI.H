/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

template<class ParcelType>
template<class TrackCloudType>
inline Foam::ThermoParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud
)
:
    ParcelType::trackingData(cloud),
    Cp_(cloud.thermo().thermo().Cp()),
    kappa_(cloud.thermo().thermo().kappa()),
    TInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.T()
        )
    ),
    CpInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            Cp_
        )
    ),
    kappaInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            kappa_
        )
    ),
    GInterp_(nullptr),
    Tc_(Zero),
    Cpc_(Zero)
{
    if (cloud.radiation())
    {
        GInterp_.reset
        (
            interpolation<scalar>::New
            (
                cloud.solution().interpolationSchemes(),
                cloud.mesh().objectRegistry::template
                    lookupObject<volScalarField>("G")
            ).ptr()
        );
    }
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::ThermoParcel<ParcelType>::trackingData::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::ThermoParcel<ParcelType>::trackingData::kappa() const
{
    return kappa_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackingData::TInterp() const
{
    return TInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackingData::CpInterp() const
{
    return CpInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackingData::kappaInterp() const
{
    return kappaInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::trackingData::GInterp() const
{
    if (!GInterp_.valid())
    {
        FatalErrorInFunction
            << "Radiation G interpolation object not set"
            << abort(FatalError);
    }

    return GInterp_();
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::trackingData::Tc() const
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::trackingData::Tc()
{
    return Tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::trackingData::Cpc() const
{
    return Cpc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::trackingData::Cpc()
{
    return Cpc_;
}


// ************************************************************************* //
