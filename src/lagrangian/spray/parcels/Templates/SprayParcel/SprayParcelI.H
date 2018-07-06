/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    sigma0_(this->dict_, 0.0),
    mu0_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    sigma0_(cp.sigma0_),
    mu0_(cp.mu0_)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    sigma0_(this->dict_, "sigma0"),
    mu0_(this->dict_, "mu0")
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar rho0,
    const scalar minParcelMass,
    const scalar youngsModulus,
    const scalar poissonsRatio,
    const scalar T0,
    const scalar TMin,
    const scalar TMax,
    const scalar Cp0,
    const scalar epsilon0,
    const scalar f0,
    const scalar Pr,
    const scalar pMin,
    const Switch& constantVolume,
    const scalar sigma0,
    const scalar mu0
)
:
    ParcelType::constantProperties
    (
        parcelTypeId,
        rhoMin,
        rho0,
        minParcelMass,
        youngsModulus,
        poissonsRatio,
        T0,
        TMin,
        TMax,
        Cp0,
        epsilon0,
        f0,
        Pr,
        pMin,
        constantVolume
    ),
    sigma0_(this->dict_, sigma0),
    mu0_(this->dict_, mu0)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    d0_(this->d()),
    position0_(this->position()),
    sigma_(0.0),
    mu_(0.0),
    liquidCore_(0.0),
    KHindex_(0.0),
    y_(0.0),
    yDot_(0.0),
    tc_(0.0),
    ms_(0.0),
    injector_(1.0),
    tMom_(great),
    user_(0.0)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    d0_(this->d()),
    position0_(this->position()),
    sigma_(0.0),
    mu_(0.0),
    liquidCore_(0.0),
    KHindex_(0.0),
    y_(0.0),
    yDot_(0.0),
    tc_(0.0),
    ms_(0.0),
    injector_(1.0),
    tMom_(great),
    user_(0.0)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const scalar liquidCore,
    const scalar KHindex,
    const scalar y,
    const scalar yDot,
    const scalar tc,
    const scalar ms,
    const scalar injector,
    const scalar tMom,
    const scalar user,
    const typename ParcelType::constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    ),
    d0_(d0),
    position0_(this->position()),
    sigma_(constProps.sigma0()),
    mu_(constProps.mu0()),
    liquidCore_(liquidCore),
    KHindex_(KHindex),
    y_(y),
    yDot_(yDot),
    tc_(tc),
    ms_(ms),
    injector_(injector),
    tMom_(tMom),
    user_(user)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::SprayParcel<ParcelType>::constantProperties::sigma0() const
{
    return sigma0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::SprayParcel<ParcelType>::constantProperties::mu0() const
{
    return mu0_.value();
}


// * * * * * * * * * * SprayParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::d0() const
{
    return d0_;
}


template<class ParcelType>
inline const Foam::vector& Foam::SprayParcel<ParcelType>::position0() const
{
    return position0_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::sigma() const
{
    return sigma_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::mu() const
{
    return mu_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::liquidCore() const
{
    return liquidCore_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::KHindex() const
{
    return KHindex_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::y() const
{
    return y_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::yDot() const
{
    return yDot_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::tc() const
{
    return tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::ms() const
{
    return ms_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::injector() const
{
    return injector_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::tMom() const
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::user() const
{
    return user_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::d0()
{
    return d0_;
}


template<class ParcelType>
inline Foam::vector& Foam::SprayParcel<ParcelType>::position0()
{
    return position0_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::sigma()
{
    return sigma_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::mu()
{
    return mu_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::liquidCore()
{
    return liquidCore_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::KHindex()
{
    return KHindex_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::y()
{
    return y_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::yDot()
{
    return yDot_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::tc()
{
    return tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::ms()
{
    return ms_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::injector()
{
    return injector_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::tMom()
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::user()
{
    return user_;
}


// ************************************************************************* //
