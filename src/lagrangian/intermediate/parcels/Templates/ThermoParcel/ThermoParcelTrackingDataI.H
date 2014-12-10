/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
template<class CloudType>
inline Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part
)
:
    ParcelType::template TrackingData<CloudType>(cloud, part),
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
    GInterp_(NULL)
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
template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::Cp() const
{
    return Cp_;
}


template<class ParcelType>
template<class CloudType>
inline const Foam::volScalarField&
Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::kappa() const
{
    return kappa_;
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::TInterp() const
{
    return TInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::CpInterp() const
{
    return CpInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::kappaInterp() const
{
    return kappaInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::GInterp() const
{
    if (!GInterp_.valid())
    {
        FatalErrorIn
        (
            "inline const Foam::interpolation<Foam::scalar>&"
            "Foam::ThermoParcel<ParcelType>::TrackingData<CloudType>::"
            "GInterp() const"
        )   << "Radiation G interpolation object not set"
            << abort(FatalError);
    }

    return GInterp_();
}


// ************************************************************************* //
