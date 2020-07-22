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
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    TDevol_(this->dict_, 0.0),
    LDevol_(this->dict_, 0.0),
    hRetentionCoeff_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    TDevol_(cp.TDevol_),
    LDevol_(cp.LDevol_),
    hRetentionCoeff_(cp.hRetentionCoeff_)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    TDevol_(this->dict_, "TDevol"),
    LDevol_(this->dict_, "LDevol"),
    hRetentionCoeff_(this->dict_, "hRetentionCoeff")
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
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
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const constantProperties& constProps
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
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0),
    canCombust_(0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::TDevol() const
{
    return TDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::LDevol() const
{
    return LDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseParcel<ParcelType>::constantProperties::
hRetentionCoeff() const
{
    scalar value = hRetentionCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorInFunction
            << "hRetentionCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }

    return value;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::
YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label
Foam::ReactingMultiphaseParcel<ParcelType>::canCombust() const
{
    return canCombust_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingMultiphaseParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label& Foam::ReactingMultiphaseParcel<ParcelType>::canCombust()
{
    return canCombust_;
}


// ************************************************************************* //
