/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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


template<class CloudType>
void Foam::MPPICCloud<CloudType>::cloudReset(MPPICCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    packingModel_.reset(c.packingModel_.ptr());
    dampingModel_.reset(c.dampingModel_.ptr());
    isotropyModel_.reset(c.isotropyModel_.ptr());
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
    cloudReset(cloudCopyPtr_());

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
    // Assign parcel ID-s
    label i = 0;
    forAllIter(typename MPPICCloud<CloudType>, *this, iter)
    {
        iter().id() = labelPair(Pstream::myProcNo(), i ++);
    }

    // Create a copy of all parcels and sources to use as a predictor
    autoPtr<MPPICCloud<CloudType>> predictorCloudPtr
    (
        static_cast<MPPICCloud<CloudType>*>
        (
            clone(this->name() + "Predictor").ptr()
        )
    );
    MPPICCloud<CloudType>& predictorCloud = predictorCloudPtr();

    // Predictor move
    predictorCloud.CloudType::move(predictorCloud, td);

    // Calculate correction velocities
    const scalar trackTime = td.trackTime();
    td.updateAverages(predictorCloud);
    predictorCloud.dampingModel().cacheFields(true);
    predictorCloud.packingModel().cacheFields(true);
    vectorField UCorr(this->size(), Zero);
    List<DynamicList<vector>> UCorrProc(Pstream::nProcs());
    List<DynamicList<label>> IDProc(Pstream::nProcs());
    forAllIter(typename MPPICCloud<CloudType>, predictorCloud, iter)
    {
        const labelPair& id = iter().id();

        const vector dU =
            predictorCloud.packingModel().velocityCorrection(iter(), trackTime)
          + predictorCloud.dampingModel().velocityCorrection(iter(), trackTime);

        if (id.first() == Pstream::myProcNo())
        {
            UCorr[id.second()] = dU;
        }
        else
        {
            UCorrProc[id.first()].append(dU);
            IDProc[id.first()].append(id.second());
        }
    }
    predictorCloud.dampingModel().cacheFields(false);
    predictorCloud.packingModel().cacheFields(false);

    // Distribute the correction velocities
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    if (Pstream::parRun())
    {
        forAll(UCorrProc, proci)
        {
            if (proci == Pstream::myProcNo()) continue;

            UOPstream os(proci, pBufs);

            os  << UCorrProc[proci] << IDProc[proci];
        }

        pBufs.finishedSends();

        forAll(UCorrProc, proci)
        {
            if (proci == Pstream::myProcNo()) continue;

            UIPstream is(proci, pBufs);

            is  >> UCorrProc[proci] >> IDProc[proci];
        }
    }
    forAll(UCorrProc, proci)
    {
        if (proci == Pstream::myProcNo()) continue;

        forAll(UCorrProc[proci], i)
        {
            UCorr[IDProc[proci][i]] = UCorrProc[proci][i];
        }
    }

    // Apply the correction velocities to the parcels
    forAllIter(typename MPPICCloud<CloudType>, *this, iter)
    {
        iter().U() += UCorr[iter().id().second()];
    }

    // Corrector
    CloudType::move(cloud, td);

    // Apply isotropy model
    td.updateAverages(cloud);
    isotropyModel_->calculate();

    // Update cell occupancy
    this->updateCellOccupancy();
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
