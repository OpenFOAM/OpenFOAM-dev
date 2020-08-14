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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    T0_(this->dict_, 0.0),
    TMin_(this->dict_, 0.0),
    TMax_(this->dict_, vGreat),
    Cp0_(this->dict_, 0.0),
    epsilon0_(this->dict_, 0.0),
    f0_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    T0_(cp.T0_),
    TMin_(cp.TMin_),
    TMax_(cp.TMax_),
    Cp0_(cp.Cp0_),
    epsilon0_(cp.epsilon0_),
    f0_(cp.f0_)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    T0_(this->dict_, "T0"),
    TMin_(this->dict_, "TMin", 200.0),
    TMax_(this->dict_, "TMax", 5000.0),
    Cp0_(this->dict_, "Cp0"),
    epsilon0_(this->dict_, "epsilon0"),
    f0_(this->dict_, "f0")
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    T_(0.0),
    Cp_(0.0)
{}


template<class ParcelType>
inline Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    T_(0.0),
    Cp_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::T0() const
{
    return T0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMin() const
{
    return TMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::TMax() const
{
    return TMax_.value();
}


template<class ParcelType>
inline void
Foam::ThermoParcel<ParcelType>::constantProperties::setTMax(const scalar TMax)
{
    TMax_.setValue(TMax);
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::Cp0() const
{
    return Cp0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::epsilon0() const
{
    return epsilon0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ThermoParcel<ParcelType>::constantProperties::f0() const
{
    return f0_.value();
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::Cp() const
{
    return Cp_;
}


template<class ParcelType>
inline Foam::scalar Foam::ThermoParcel<ParcelType>::hs() const
{
    return Cp_*(T_ - 298.15);
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::T()
{
    return T_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ThermoParcel<ParcelType>::Cp()
{
    return Cp_;
}


// ************************************************************************* //
