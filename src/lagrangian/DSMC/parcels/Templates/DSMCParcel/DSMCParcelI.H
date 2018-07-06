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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::constantProperties::constantProperties()
:
    mass_(0),
    d_(0)
{}


template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(readScalar(dict.lookup("mass"))),
    d_(readScalar(dict.lookup("diameter"))),
    internalDegreesOfFreedom_
    (
        readInt(dict.lookup("internalDegreesOfFreedom"))
    ),
    omega_(readScalar(dict.lookup("omega")))
{}


template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::DSMCParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const vector& U,
    const scalar Ei,
    const label typeId
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    U_(U),
    Ei_(Ei),
    typeId_(typeId)
{}


template<class ParcelType>
inline Foam::DSMCParcel<ParcelType>::DSMCParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const vector& U,
    const scalar Ei,
    const label typeId
)
:
    ParcelType(mesh, position, celli),
    U_(U),
    Ei_(Ei),
    typeId_(typeId)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::DSMCParcel<ParcelType>::constantProperties::mass() const
{
    return mass_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCParcel<ParcelType>::constantProperties::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCParcel<ParcelType>::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}


template<class ParcelType>
inline Foam::direction
Foam::DSMCParcel<ParcelType>::constantProperties::internalDegreesOfFreedom()
const
{
    return internalDegreesOfFreedom_;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCParcel<ParcelType>::constantProperties::omega() const
{
    return omega_;
}


// * * * * * * * * * * DSMCParcel Member Functions  * * * * * * * * * * //

template<class ParcelType>
inline Foam::label Foam::DSMCParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline const Foam::vector& Foam::DSMCParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCParcel<ParcelType>::Ei() const
{
    return Ei_;
}


template<class ParcelType>
inline Foam::vector& Foam::DSMCParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar& Foam::DSMCParcel<ParcelType>::Ei()
{
    return Ei_;
}


// ************************************************************************* //
