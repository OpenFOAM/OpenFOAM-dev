/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "CollidingCloud.H"
#include "subCycleTime.H"

#include "CollisionModel.H"
#include "NoCollision.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::CollidingCloud<CloudType>::setModels()
{
    collisionModel_.reset
    (
        CollisionModel<CollidingCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::cloudReset(CollidingCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    collisionModel_.reset(c.collisionModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    const bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    constProps_(this->particleProperties()),
    collisionModel_(nullptr)
{
    setModels();

    if (readFields)
    {
        parcelType::readFields(*this);
        this->deleteLostParticles();
    }

    if
    (
        this->solution().steadyState()
    && !isType<NoCollision<CollidingCloud<CloudType>>>(collisionModel_())
    )
    {
        FatalErrorInFunction
            << "Collision modelling not currently available "
            << "for steady state calculations" << exit(FatalError);
    }
}


template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo,
    const bool readFields
)
:
    CollidingCloud(cloudName, rho, U, carrierThermo.mu(), g, readFields)
{}


template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    CollidingCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    collisionModel_(c.collisionModel_->clone())
{}


template<class CloudType>
Foam::CollidingCloud<CloudType>::CollidingCloud
(
    const fvMesh& mesh,
    const word& name,
    const CollidingCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    collisionModel_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CollidingCloud<CloudType>::~CollidingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CollidingCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<CollidingCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
template<class TrackCloudType>
void  Foam::CollidingCloud<CloudType>::motion
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    // Sympletic leapfrog integration of particle forces:
    // + apply half deltaV with stored force
    // + move positions with new velocity
    // + calculate forces in new position
    // + apply half deltaV with new force

    const label nSubCycles = collision().nSubCycles();

    if (nSubCycles > 1)
    {
        Info<< "    " << nSubCycles << " move-collide subCycles" << endl;
    }

    for (label subCyclei = 0; subCyclei < nSubCycles; ++ subCyclei)
    {
        td.stepFractionRange() =
            Pair<scalar>
            (
                scalar(subCyclei)/nSubCycles,
                scalar(subCyclei + 1)/nSubCycles
            );

        td.part() = parcelType::trackingData::tpVelocityHalfStep;
        CloudType::move(cloud, td);

        td.part() = parcelType::trackingData::tpLinearTrack;
        CloudType::move(cloud, td);

        // td.part() = parcelType::trackingData::tpRotationalTrack;
        // CloudType::move(cloud, td);

        this->updateCellOccupancy();

        this->collision().collide();

        td.part() = parcelType::trackingData::tpVelocityHalfStep;
        CloudType::move(cloud, td);
    }

    td.stepFractionRange() = Pair<scalar>(0, 1);
}


template<class CloudType>
void Foam::CollidingCloud<CloudType>::info()
{
    CloudType::info();

    scalar rotationalKineticEnergy = rotationalKineticEnergyOfSystem();
    reduce(rotationalKineticEnergy, sumOp<scalar>());

    Info<< "    Rotational kinetic energy       = "
        << rotationalKineticEnergy << nl;
}


// ************************************************************************* //
