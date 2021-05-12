/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "MPPICCloud.H"
#include "NoPacking.H"
#include "ParticleStressModel.H"
#include "NoDamping.H"
#include "NoIsotropy.H"
#include "TimeScaleModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::MPPICCloud<CloudType>::setModels()
{
    packingModel_.reset
    (
        PackingModel<MPPICCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    dampingModel_.reset
    (
        DampingModel<MPPICCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
    isotropyModel_.reset
    (
        IsotropyModel<MPPICCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MPPICCloud<CloudType>::MPPICCloud
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
    packingModel_(nullptr),
    dampingModel_(nullptr),
    isotropyModel_(nullptr)
{
    if (this->solution().steadyState())
    {
        FatalErrorInFunction
            << "MPPIC modelling not available for steady state calculations"
            << exit(FatalError);
    }

    setModels();

    if (readFields)
    {
        parcelType::readFields(*this);
        this->deleteLostParticles();
    }
}


template<class CloudType>
Foam::MPPICCloud<CloudType>::MPPICCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo,
    const bool readFields
)
:
    MPPICCloud(cloudName, rho, U, carrierThermo.mu(), g, readFields)
{}


template<class CloudType>
Foam::MPPICCloud<CloudType>::MPPICCloud
(
    MPPICCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    packingModel_(c.packingModel_->clone()),
    dampingModel_(c.dampingModel_->clone()),
    isotropyModel_(c.isotropyModel_->clone())
{}


template<class CloudType>
Foam::MPPICCloud<CloudType>::MPPICCloud
(
    const fvMesh& mesh,
    const word& name,
    const MPPICCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    packingModel_(nullptr),
    dampingModel_(nullptr),
    isotropyModel_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MPPICCloud<CloudType>::~MPPICCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MPPICCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<MPPICCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::MPPICCloud<CloudType>::restoreState()
{
    this->cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::MPPICCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::MPPICCloud<CloudType>::motion
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    // Momentum
    // ~~~~~~~~

    // force calculation and tracking
    td.part() = parcelType::trackingData::tpPredictTrack;
    CloudType::move(cloud, td, this->db().time().deltaTValue());


    // Preliminary
    // ~~~~~~~~~~~

    // switch forces off so they are not applied in corrector steps
    this->forces().setCalcNonCoupled(false);
    this->forces().setCalcCoupled(false);


    // Damping
    // ~~~~~~~

    if (!isType<DampingModels::NoDamping<CloudType>>(dampingModel_()))
    {
        if (this->mesh().moving())
        {
            FatalErrorInFunction
                << "MPPIC damping modelling does not support moving meshes."
                << exit(FatalError);
        }

        // update averages
        td.updateAverages(cloud);

        // memory allocation and eulerian calculations
        dampingModel_->cacheFields(true);

        // calculate the damping velocity corrections without moving the parcels
        td.part() = parcelType::trackingData::tpDampingNoTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());

        // correct the parcel positions and velocities
        td.part() = parcelType::trackingData::tpCorrectTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());

        // finalise and free memory
        dampingModel_->cacheFields(false);
    }


    // Packing
    // ~~~~~~~

    if (!isType<PackingModels::NoPacking<CloudType>>(packingModel_()))
    {
        if (this->mesh().moving())
        {
            FatalErrorInFunction
                << "MPPIC packing modelling does not support moving meshes."
                << exit(FatalError);
        }

        // same procedure as for damping
        td.updateAverages(cloud);
        packingModel_->cacheFields(true);
        td.part() = parcelType::trackingData::tpPackingNoTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());
        td.part() = parcelType::trackingData::tpCorrectTrack;
        CloudType::move(cloud, td, this->db().time().deltaTValue());
        packingModel_->cacheFields(false);
    }


    // Isotropy
    // ~~~~~~~~

    if (!isType<IsotropyModels::NoIsotropy<CloudType>>(isotropyModel_()))
    {
        // update averages
        td.updateAverages(cloud);

        // apply isotropy model
        isotropyModel_->calculate();
    }


    // Final
    // ~~~~~

    // update cell occupancy
    this->updateCellOccupancy();

    // switch forces back on
    this->forces().setCalcNonCoupled(true);
    this->forces().setCalcCoupled(this->solution().coupled());
}


template<class CloudType>
void Foam::MPPICCloud<CloudType>::info()
{
    CloudType::info();

    tmp<volScalarField> alpha = this->theta();

    const scalar alphaMin = gMin(alpha().primitiveField());
    const scalar alphaMax = gMax(alpha().primitiveField());

    Info<< "    Min cell volume fraction        = " << alphaMin << endl;
    Info<< "    Max cell volume fraction        = " << alphaMax << endl;

    if (alphaMax < small)
    {
        return;
    }

    scalar nMin = great;

    forAll(this->mesh().cells(), celli)
    {
        const label n = this->cellOccupancy()[celli].size();

        if (n > 0)
        {
            const scalar nPack = n*alphaMax/alpha()[celli];

            if (nPack < nMin)
            {
                nMin = nPack;
            }
        }
    }

    reduce(nMin, minOp<scalar>());

    Info<< "    Min dense number of parcels     = " << nMin << endl;
}


// ************************************************************************* //
