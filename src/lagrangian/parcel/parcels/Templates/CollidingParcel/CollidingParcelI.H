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
inline Foam::CollidingParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    youngsModulus_(this->dict_, 0.0),
    poissonsRatio_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::CollidingParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    youngsModulus_(cp.youngsModulus_),
    poissonsRatio_(cp.poissonsRatio_)
{}


template<class ParcelType>
inline Foam::CollidingParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    youngsModulus_(this->dict_, "youngsModulus"),
    poissonsRatio_(this->dict_, "poissonsRatio")
{}


template<class ParcelType>
inline Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    f_(Zero),
    angularMomentum_(Zero),
    torque_(Zero),
    collisionRecords_()
{}


template<class ParcelType>
inline Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli
)
:
    ParcelType(owner, position, celli),
    f_(Zero),
    angularMomentum_(Zero),
    torque_(Zero),
    collisionRecords_()
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::CollidingParcel<ParcelType>::constantProperties::youngsModulus() const
{
    return youngsModulus_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::CollidingParcel<ParcelType>::constantProperties::poissonsRatio() const
{
    return poissonsRatio_.value();
}


// * * * * * * * * * * CollidingParcel Member Functions  * * * * * * * * * * //

template<class ParcelType>
inline const Foam::vector& Foam::CollidingParcel<ParcelType>::f() const
{
    return f_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::CollidingParcel<ParcelType>::angularMomentum() const
{
    return angularMomentum_;
}


template<class ParcelType>
inline const Foam::vector& Foam::CollidingParcel<ParcelType>::torque() const
{
    return torque_;
}


template<class ParcelType>
inline const Foam::collisionRecordList&
Foam::CollidingParcel<ParcelType>::collisionRecords() const
{
    return collisionRecords_;
}


template<class ParcelType>
inline Foam::vector& Foam::CollidingParcel<ParcelType>::f()
{
    return f_;
}


template<class ParcelType>
inline Foam::vector& Foam::CollidingParcel<ParcelType>::angularMomentum()
{
    return angularMomentum_;
}


template<class ParcelType>
inline Foam::vector& Foam::CollidingParcel<ParcelType>::torque()
{
    return torque_;
}


template<class ParcelType>
inline Foam::collisionRecordList&
Foam::CollidingParcel<ParcelType>::collisionRecords()
{
    return collisionRecords_;
}


template<class ParcelType>
inline Foam::vector Foam::CollidingParcel<ParcelType>::omega() const
{
    return angularMomentum_/this->momentOfInertia();
}


// ************************************************************************* //
